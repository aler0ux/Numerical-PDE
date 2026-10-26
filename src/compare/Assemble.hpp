#ifndef PDE_COMPARE_ASSEMBLE_HPP
#define PDE_COMPARE_ASSEMBLE_HPP

// Two assemblers for the same problem, -Laplacian(u) = f with u = 0 on the
// boundary, on the same Grid2d:
//
//   assembleFDM  -- standard 5-point finite-difference stencil.
//   assembleFEM  -- linear (P1) finite elements on a triangulation of the
//                   grid's active cells.
//
// Both build a plain list of (row, col, value) triplets plus a load
// vector, then go through the *same* applyDirichletZero() pass before
// becoming a CSRMatrix. That shared pass is what keeps both systems
// symmetric positive definite (so the one Conjugate-Gradient solver in
// ConjugateGradient.hpp serves both), and it is also the only place either
// assembler has to think about the boundary at all -- interior assembly
// can stay oblivious to it.

#include <cmath>
#include <vector>

#include "CSRMatrix.hpp"
#include "Grid2d.hpp"
#include "SourceFunctions.hpp"

namespace cmp {

struct System {
    CSRMatrix A;
    std::vector<double> b;
};

// Zero out every row and column belonging to a boundary DOF, replace the
// row with an identity entry, and zero the corresponding right-hand side.
// Because the boundary condition is homogeneous (u = 0), there is no
// correction term to move to the right-hand side of the other rows --
// dropping those columns *is* the correction.
inline void applyDirichletZero(std::vector<Triplet>& triplets, std::vector<double>& b, const Grid2d& grid) {
    std::vector<char> isBoundaryDof(static_cast<size_t>(grid.numDofs()), 0);
    for (int i = 0; i <= grid.n(); ++i)
        for (int j = 0; j <= grid.n(); ++j)
            if (grid.active(i, j) && grid.boundary(i, j)) isBoundaryDof[static_cast<size_t>(grid.dof(i, j))] = 1;

    std::vector<Triplet> kept;
    kept.reserve(triplets.size());
    for (const Triplet& t : triplets) {
        if (isBoundaryDof[static_cast<size_t>(t.row)] || isBoundaryDof[static_cast<size_t>(t.col)]) continue;
        kept.push_back(t);
    }
    for (int d = 0; d < grid.numDofs(); ++d) {
        if (!isBoundaryDof[static_cast<size_t>(d)]) continue;
        kept.push_back(Triplet{d, d, 1.0});
        b[static_cast<size_t>(d)] = 0.0;
    }
    triplets.swap(kept);
}

// ---------------------------------------------------------------------------
// Finite differences: standard 5-point Laplacian, h^2 uniform spacing.
// ---------------------------------------------------------------------------
inline System assembleFDM(const Grid2d& grid, const pde::Expr& f) {
    const int n = grid.numDofs();
    const double h2 = grid.h() * grid.h();
    std::vector<Triplet> triplets;
    std::vector<double> b(static_cast<size_t>(n), 0.0);

    for (int i = 0; i <= grid.n(); ++i) {
        for (int j = 0; j <= grid.n(); ++j) {
            if (!grid.active(i, j)) continue;
            const int d = grid.dof(i, j);
            b[static_cast<size_t>(d)] = evalF(f, grid.x(i), grid.y(j));

            if (grid.boundary(i, j)) continue;  // row replaced by applyDirichletZero()

            // Non-boundary active nodes are guaranteed (by Grid2d::build())
            // to have all four neighbours active, so this never reaches
            // outside the domain.
            triplets.push_back(Triplet{d, d, 4.0 / h2});
            triplets.push_back(Triplet{d, grid.dof(i - 1, j), -1.0 / h2});
            triplets.push_back(Triplet{d, grid.dof(i + 1, j), -1.0 / h2});
            triplets.push_back(Triplet{d, grid.dof(i, j - 1), -1.0 / h2});
            triplets.push_back(Triplet{d, grid.dof(i, j + 1), -1.0 / h2});
        }
    }

    applyDirichletZero(triplets, b, grid);
    return System{CSRMatrix::fromTriplets(n, std::move(triplets)), std::move(b)};
}

// ---------------------------------------------------------------------------
// Finite elements: linear (P1) triangles over the grid's active cells.
// ---------------------------------------------------------------------------
namespace detail {

struct Vertex {
    int dof;
    double x, y;
};

// One triangle's contribution to the stiffness matrix and load vector.
// Stiffness: K_ij = integral(grad(lambda_i) . grad(lambda_j)) dA, which for
// a linear triangle has the closed form (b_i b_j + c_i c_j) / (4 * Area).
// Load: single centroid-point quadrature, f(centroid) * Area / 3 per
// vertex -- Area/3 is the exact integral of a single linear basis function
// over the triangle, so this is exact for constant f and a simple, standard
// one-point approximation otherwise.
inline void addTriangle(const Vertex& v0, const Vertex& v1, const Vertex& v2, const pde::Expr& f,
                         std::vector<Triplet>& triplets, std::vector<double>& b) {
    const Vertex verts[3] = {v0, v1, v2};
    const double area2 = (v1.x - v0.x) * (v2.y - v0.y) - (v2.x - v0.x) * (v1.y - v0.y);
    const double area = 0.5 * std::fabs(area2);
    if (area < 1e-14) return;  // degenerate triangle, nothing to add

    // b_i, c_i are the (unnormalised) gradient components of barycentric
    // coordinate lambda_i: b_i = y_{i+1} - y_{i-1}, c_i = x_{i-1} - x_{i+1}
    // (indices mod 3).
    double bc[3], cc[3];
    for (int i = 0; i < 3; ++i) {
        const Vertex& vNext = verts[(i + 1) % 3];
        const Vertex& vPrev = verts[(i + 2) % 3];
        bc[i] = vNext.y - vPrev.y;
        cc[i] = vPrev.x - vNext.x;
    }

    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            const double k = (bc[i] * bc[j] + cc[i] * cc[j]) / (4.0 * area);
            triplets.push_back(Triplet{verts[i].dof, verts[j].dof, k});
        }
    }

    const double cx = (v0.x + v1.x + v2.x) / 3.0;
    const double cy = (v0.y + v1.y + v2.y) / 3.0;
    const double fc = evalF(f, cx, cy);
    for (int i = 0; i < 3; ++i) b[static_cast<size_t>(verts[i].dof)] += fc * area / 3.0;
}

}  // namespace detail

inline System assembleFEM(const Grid2d& grid, const pde::Expr& f) {
    const int n = grid.numDofs();
    std::vector<Triplet> triplets;
    std::vector<double> b(static_cast<size_t>(n), 0.0);

    auto vertex = [&](int i, int j) -> detail::Vertex { return {grid.dof(i, j), grid.x(i), grid.y(j)}; };

    // Triangulate every grid cell from its active corners. A cell with all
    // four corners active splits along the BL-TR diagonal into two
    // triangles, the standard structured-mesh triangulation. A cell with
    // exactly three active corners still contributes the one triangle
    // those three corners form, so the FEM mesh's footprint tracks the
    // domain boundary a little more closely than skipping partial cells
    // would. Fewer than three active corners contributes nothing.
    for (int i = 0; i < grid.n(); ++i) {
        for (int j = 0; j < grid.n(); ++j) {
            const bool bl = grid.active(i, j);
            const bool br = grid.active(i + 1, j);
            const bool tr = grid.active(i + 1, j + 1);
            const bool tl = grid.active(i, j + 1);
            const int activeCount = (bl ? 1 : 0) + (br ? 1 : 0) + (tr ? 1 : 0) + (tl ? 1 : 0);
            if (activeCount < 3) continue;

            if (activeCount == 4) {
                detail::addTriangle(vertex(i, j), vertex(i + 1, j), vertex(i + 1, j + 1), f, triplets, b);
                detail::addTriangle(vertex(i, j), vertex(i + 1, j + 1), vertex(i, j + 1), f, triplets, b);
            } else {
                std::vector<detail::Vertex> corners;
                if (bl) corners.push_back(vertex(i, j));
                if (br) corners.push_back(vertex(i + 1, j));
                if (tr) corners.push_back(vertex(i + 1, j + 1));
                if (tl) corners.push_back(vertex(i, j + 1));
                detail::addTriangle(corners[0], corners[1], corners[2], f, triplets, b);
            }
        }
    }

    // Any active DOF that no triangle ever touched (an isolated node, e.g.
    // a single active corner diagonally adjacent to the rest of the
    // domain) would otherwise leave a structurally zero row. Pin it to
    // zero exactly like a boundary node; it is geometrically on the edge
    // of what the FEM mesh can represent either way.
    std::vector<char> touched(static_cast<size_t>(n), 0);
    for (const Triplet& t : triplets) touched[static_cast<size_t>(t.row)] = 1;
    for (int i = 0; i <= grid.n(); ++i) {
        for (int j = 0; j <= grid.n(); ++j) {
            if (!grid.active(i, j)) continue;
            const int d = grid.dof(i, j);
            if (!touched[static_cast<size_t>(d)]) triplets.push_back(Triplet{d, d, 1.0});
        }
    }

    applyDirichletZero(triplets, b, grid);
    return System{CSRMatrix::fromTriplets(n, std::move(triplets)), std::move(b)};
}

}  // namespace cmp

#endif  // PDE_COMPARE_ASSEMBLE_HPP

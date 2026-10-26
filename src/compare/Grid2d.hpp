#ifndef PDE_COMPARE_GRID2D_HPP
#define PDE_COMPARE_GRID2D_HPP

// A structured (N+1) x (N+1) grid over [0,1] x [0,1], classified against
// one of the Domains.hpp presets (each a pde::Domain2d CSG shape -- see
// that file). This single grid is shared by both
// discretizations in this module -- the FDM stencil and the FEM mesh are
// both built on exactly the same node set and the same degree-of-freedom
// numbering, which is what makes comparing the two solutions node-for-node
// meaningful instead of needing to interpolate between two different
// meshes.
//
// Boundary nodes keep a DOF (same choice the roadmap makes for the general
// assembler): their row becomes "u = 0" instead of being eliminated, which
// keeps both assemblers simple and uniform.

#include <vector>

#include "Domains.hpp"

namespace cmp {

class Grid2d {
public:
    Grid2d(DomainKind kind, int n)
        : kind_(kind), n_(n), h_(1.0 / static_cast<double>(n)), domain_(buildDomain2d(kind)) {
        const int side = n_ + 1;
        active_.assign(static_cast<size_t>(side) * static_cast<size_t>(side), false);
        boundary_.assign(static_cast<size_t>(side) * static_cast<size_t>(side), false);
        dof_.assign(static_cast<size_t>(side) * static_cast<size_t>(side), -1);
        build();
    }

    int n() const { return n_; }          // cells per axis; side() = n()+1 nodes per axis
    int side() const { return n_ + 1; }
    double h() const { return h_; }
    DomainKind kind() const { return kind_; }

    double x(int i) const { return static_cast<double>(i) * h_; }
    double y(int j) const { return static_cast<double>(j) * h_; }

    int index(int i, int j) const { return i * side() + j; }

    bool active(int i, int j) const { return active_[static_cast<size_t>(index(i, j))]; }
    bool boundary(int i, int j) const { return boundary_[static_cast<size_t>(index(i, j))]; }
    int dof(int i, int j) const { return dof_[static_cast<size_t>(index(i, j))]; }

    int numDofs() const { return numDofs_; }

private:
    DomainKind kind_;
    int n_;
    double h_;
    pde::Domain2d domain_;
    std::vector<char> active_;
    std::vector<char> boundary_;
    std::vector<int> dof_;
    int numDofs_ = 0;

    bool inBounds(int i, int j) const { return i >= 0 && i <= n_ && j >= 0 && j <= n_; }

    bool activeAt(int i, int j) const {
        if (!inBounds(i, j)) return false;
        return domain_.inDomain(x(i), y(j));
    }

    void build() {
        const int side_ = side();
        for (int i = 0; i <= n_; ++i) {
            for (int j = 0; j <= n_; ++j) {
                active_[static_cast<size_t>(index(i, j))] = activeAt(i, j);
            }
        }
        // A node is a boundary node if it is active and sits on the edge of
        // the bounding grid, or has a neighbour that is not active -- i.e.
        // it is where the discrete domain meets its own edge.
        for (int i = 0; i <= n_; ++i) {
            for (int j = 0; j <= n_; ++j) {
                if (!active(i, j)) continue;
                const bool onGridEdge = (i == 0 || i == n_ || j == 0 || j == n_);
                const bool neighbourMissing = !activeAt(i - 1, j) || !activeAt(i + 1, j) ||
                                               !activeAt(i, j - 1) || !activeAt(i, j + 1);
                boundary_[static_cast<size_t>(index(i, j))] = onGridEdge || neighbourMissing;
            }
        }
        // Number every active node (interior and boundary alike).
        int next = 0;
        for (int i = 0; i <= n_; ++i) {
            for (int j = 0; j <= n_; ++j) {
                if (active(i, j)) dof_[static_cast<size_t>(index(i, j))] = next++;
            }
        }
        numDofs_ = next;
        (void)side_;
    }
};

}  // namespace cmp

#endif  // PDE_COMPARE_GRID2D_HPP

#ifndef PDE_COMPARE_CONJUGATE_GRADIENT_HPP
#define PDE_COMPARE_CONJUGATE_GRADIENT_HPP

// The "simple iterative solver" requested for this project: preconditioned
// Conjugate Gradient on a CSR matrix.
//
// CG applies because every system built in this module (FDM 5-point
// Laplacian, FEM P1 stiffness) is symmetric positive definite once the
// zero-Dirichlet boundary is eliminated the way Assemble.hpp does it: the
// boundary rows/columns become an identity block, and the remaining
// interior block is the usual SPD Poisson operator. No general-purpose
// solver (LU, BiCGSTAB, ...) is needed for this comparison.
//
// Preconditioning is Jacobi (diagonal) -- the cheapest preconditioner that
// still noticeably cuts the iteration count, and it needs nothing beyond
// the diagonal CSRMatrix already exposes.

#include <cmath>
#include <vector>

#include "CSRMatrix.hpp"

namespace cmp {

struct CGResult {
    std::vector<double> x;
    std::vector<double> residualHistory;  // ||r||_2 / ||b||_2 after each iteration
    int iterations = 0;
    bool converged = false;
};

inline double dot(const std::vector<double>& a, const std::vector<double>& b) {
    double s = 0.0;
    for (size_t i = 0; i < a.size(); ++i) s += a[i] * b[i];
    return s;
}

inline double norm2(const std::vector<double>& a) { return std::sqrt(dot(a, a)); }

// Solves A x = b. x starts at the zero vector, which is the natural guess
// for these zero-Dirichlet problems.
inline CGResult conjugateGradient(const CSRMatrix& A, const std::vector<double>& b,
                                   double tol = 1e-10, int maxIterations = 20000) {
    const int n = A.rows();
    CGResult result;
    result.x.assign(static_cast<size_t>(n), 0.0);

    const std::vector<double> diag = A.diagonal();
    auto applyPreconditioner = [&](const std::vector<double>& r, std::vector<double>& z) {
        z.resize(r.size());
        for (int i = 0; i < n; ++i) {
            const double d = diag[static_cast<size_t>(i)];
            z[static_cast<size_t>(i)] = (d != 0.0) ? r[static_cast<size_t>(i)] / d : r[static_cast<size_t>(i)];
        }
    };

    const double bNorm = std::max(norm2(b), 1e-300);

    std::vector<double> r = b;  // r0 = b - A*x0 = b, since x0 = 0
    std::vector<double> z;
    applyPreconditioner(r, z);
    std::vector<double> p = z;
    double rz = dot(r, z);

    std::vector<double> Ap(static_cast<size_t>(n));

    for (int iter = 0; iter < maxIterations; ++iter) {
        A.matvec(p, Ap);
        const double pAp = dot(p, Ap);
        const double alpha = (pAp != 0.0) ? rz / pAp : 0.0;

        for (int i = 0; i < n; ++i) {
            result.x[static_cast<size_t>(i)] += alpha * p[static_cast<size_t>(i)];
            r[static_cast<size_t>(i)] -= alpha * Ap[static_cast<size_t>(i)];
        }

        const double relResidual = norm2(r) / bNorm;
        result.residualHistory.push_back(relResidual);
        result.iterations = iter + 1;

        if (relResidual < tol) {
            result.converged = true;
            break;
        }

        applyPreconditioner(r, z);
        const double rzNew = dot(r, z);
        const double beta = (rz != 0.0) ? rzNew / rz : 0.0;
        for (int i = 0; i < n; ++i) p[static_cast<size_t>(i)] = z[static_cast<size_t>(i)] + beta * p[static_cast<size_t>(i)];
        rz = rzNew;
    }

    return result;
}

}  // namespace cmp

#endif  // PDE_COMPARE_CONJUGATE_GRADIENT_HPP

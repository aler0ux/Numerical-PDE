#ifndef PDE_CORE_OPERATOR_HPP
#define PDE_CORE_OPERATOR_HPP

// The single canonical form for every linear PDE up to second order.
//
//     m2(x) d2u/dt2  +  m1(x) du/dt   =   L[u]  +  f(x,t)
//
// where L is a sum of terms, each a coefficient times one derivative of u:
//
//     L[u] = sum over alpha of  c_alpha(x,t) * d^alpha u,     order(alpha) <= 2
//
// Everything the solver supports is a choice of coefficients:
//
//     Laplace / Poisson     m2=m1=0,  c_{d2/dxi2}=1
//     Helmholtz             m2=m1=0,  c_{d2/dxi2}=1,  c_u=k^2
//     Heat / diffusion      m1=1,     c_{d2/dxi2}=kappa
//     Advection-diffusion   m1=1,     c_{d2/dxi2}=kappa,  c_{d/dxi}=-v_i
//     Wave                  m2=1,     c_{d2/dxi2}=c^2
//
// This replaces the earlier GeneralizedFirstOrder and GeneralizedSecondOrder
// structs. Their named coefficients (A, B, C, D, E, F) dead-ended at two
// dimensions and could not express both a time derivative and a full set of
// second-order terms; keying by MultiIndex removes both limits. There is
// deliberately no `general2d(A, B, C, D, E, F)` builder here -- that would
// reintroduce exactly the naming problem this form exists to remove. Use add().
//
// Sign convention: time derivatives are on the left, everything else on the
// right, and coefficients are stored exactly as written above with no hidden
// factors. A steady problem is therefore 0 = L[u] + f, so `u_xx = rho` is
// stored with source -rho. PDE::poisson() handles that flip for you.

#include <algorithm>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "Expr.hpp"
#include "MultiIndex.hpp"

namespace pde {

// One coefficient-times-derivative pair.
struct Term {
    MultiIndex alpha;
    Expr       coefficient;
};

// ---------------------------------------------------------------------------
// The spatial operator L, plus its source term f.
// ---------------------------------------------------------------------------
class LinearOperator {
public:
    explicit LinearOperator(int dim) : dim_(dim) {
        if (dim < 1 || dim > kMaxDim)
            throw std::invalid_argument("pde::LinearOperator: dim must be in [1, " +
                                        std::to_string(kMaxDim) + "]");
    }

    int dim() const { return dim_; }

    // Add `coeff * d^alpha u`. Adding the same multi-index twice accumulates
    // rather than duplicating, so the term list stays canonical: at most one
    // entry per derivative, highest order first.
    LinearOperator& add(const MultiIndex& alpha, Expr coeff) {
        if (alpha.order() > 2)
            throw std::invalid_argument("pde::LinearOperator: order " +
                                        std::to_string(alpha.order()) +
                                        " exceeds the supported maximum of 2");
        for (int i = dim_; i < kMaxDim; ++i)
            if (alpha[i] != 0)
                throw std::invalid_argument("pde::LinearOperator: multi-index " +
                                            alpha.toString() + " uses an axis beyond dim=" +
                                            std::to_string(dim_));

        auto cmp = [](const Term& term, const MultiIndex& key) {
            return PrincipalFirst{}(term.alpha, key);
        };
        auto it = std::lower_bound(terms_.begin(), terms_.end(), alpha, cmp);

        if (it != terms_.end() && it->alpha == alpha) {
            it->coefficient = it->coefficient + std::move(coeff);
            if (isZero(it->coefficient)) terms_.erase(it);
        } else if (!isZero(coeff)) {
            terms_.insert(it, Term{alpha, std::move(coeff)});
        }
        return *this;
    }

    LinearOperator& setSource(Expr f) {
        source_ = std::move(f);
        return *this;
    }

    const std::vector<Term>& terms() const { return terms_; }
    const Expr&              source() const { return source_; }

    // Coefficient of d^alpha u, or nullptr if the operator has no such term.
    const Expr* coefficient(const MultiIndex& alpha) const {
        for (const Term& term : terms_)
            if (term.alpha == alpha) return &term.coefficient;
        return nullptr;
    }

    // Highest derivative order present. 2 => needs second-order stencils.
    int maxOrder() const { return terms_.empty() ? 0 : terms_.front().alpha.order(); }

    std::string toString() const {
        std::string s;
        for (const Term& term : terms_) {
            Expr       magnitude;
            const bool negative = splitSign(term.coefficient, magnitude);
            if (s.empty()) {
                if (negative) s += "-";
            } else {
                s += negative ? " - " : " + ";
            }
            // Printed in multiplicative context, so a sum such as
            // `1 + sin(x)` brackets instead of swallowing the derivative.
            if (!isOne(magnitude)) s += toStringAt(magnitude, 2) + "*";
            s += term.alpha.toString();
        }
        if (s.empty()) s = "0";
        if (!isZero(source_)) {
            Expr       magnitude;
            const bool negative = splitSign(source_, magnitude);
            // A negated sum must keep its brackets: `- (a + b)`, not `- a + b`.
            s += (negative ? " - " : " + ") + toStringAt(magnitude, negative ? 2 : 1);
        }
        return s;
    }

    // -- common operators ---------------------------------------------------

    // kappa * (d2/dx2 + d2/dy2 + ...)
    static LinearOperator laplacian(int dim, Expr kappa = Number{1.0}) {
        LinearOperator L(dim);
        for (int i = 0; i < dim; ++i) L.add(MultiIndex::d2(i), kappa);
        return L;
    }

    // laplacian + k^2 * u
    static LinearOperator helmholtz(int dim, Expr k2) {
        LinearOperator L = laplacian(dim);
        L.add(MultiIndex::u(), std::move(k2));
        return L;
    }

    // kappa * laplacian - v . grad
    static LinearOperator advectionDiffusion(int dim, Expr kappa, std::vector<Expr> velocity) {
        if (static_cast<int>(velocity.size()) != dim)
            throw std::invalid_argument("pde::LinearOperator::advectionDiffusion: velocity has " +
                                        std::to_string(velocity.size()) + " components, expected " +
                                        std::to_string(dim));
        LinearOperator L = laplacian(dim, std::move(kappa));
        for (int i = 0; i < dim; ++i) L.add(MultiIndex::d(i), -std::move(velocity[static_cast<size_t>(i)]));
        return L;
    }

private:
    int               dim_;
    std::vector<Term> terms_;
    Expr              source_ = Number{0.0};
};

// ---------------------------------------------------------------------------
// How the problem behaves in time. This drives which solver path runs, so it
// classifies the time structure only -- it makes no claim about ellipticity,
// which depends on the eigenvalues of the principal part and is pointwise when
// coefficients vary.
// ---------------------------------------------------------------------------
enum class TimeKind {
    Steady,       // boundary value problem, solved once
    FirstOrder,   // du/dt present: theta-method
    SecondOrder   // d2u/dt2 present: reduce to a first-order system
};

// ---------------------------------------------------------------------------
// A complete equation.
// ---------------------------------------------------------------------------
struct PDE {
    LinearOperator spatial;
    Expr           mass1 = Number{0.0};  // coefficient of du/dt
    Expr           mass2 = Number{0.0};  // coefficient of d2u/dt2

    explicit PDE(LinearOperator L) : spatial(std::move(L)) {}

    int dim() const { return spatial.dim(); }

    TimeKind timeKind() const {
        if (!isZero(mass2)) return TimeKind::SecondOrder;
        if (!isZero(mass1)) return TimeKind::FirstOrder;
        return TimeKind::Steady;
    }

    bool isSteady() const { return timeKind() == TimeKind::Steady; }

    std::string toString() const {
        auto timeTerm = [](const Expr& m, const char* deriv) -> std::string {
            if (isZero(m)) return {};
            return isOne(m) ? std::string(deriv) : (pde::toString(m) + "*" + deriv);
        };
        std::string lhs = timeTerm(mass2, "d2u/dt2");
        const std::string first = timeTerm(mass1, "du/dt");
        if (!first.empty()) lhs = lhs.empty() ? first : lhs + " + " + first;
        if (lhs.empty()) lhs = "0";
        return lhs + " = " + spatial.toString();
    }

    // -- common equations ---------------------------------------------------

    static PDE steady(LinearOperator L) { return PDE(std::move(L)); }

    // laplacian(u) = rho   (note the source flip described at the top)
    static PDE poisson(int dim, Expr rho) {
        LinearOperator L = LinearOperator::laplacian(dim);
        L.setSource(-std::move(rho));
        return PDE(std::move(L));
    }

    // du/dt = kappa * laplacian(u)
    static PDE heat(int dim, Expr kappa) {
        PDE p(LinearOperator::laplacian(dim, std::move(kappa)));
        p.mass1 = Number{1.0};
        return p;
    }

    // d2u/dt2 = speed^2 * laplacian(u)
    static PDE wave(int dim, Expr speed) {
        PDE p(LinearOperator::laplacian(dim, speed * speed));
        p.mass2 = Number{1.0};
        return p;
    }

    // du/dt = kappa * laplacian(u) - v . grad(u)
    static PDE advectionDiffusion(int dim, Expr kappa, std::vector<Expr> velocity) {
        PDE p(LinearOperator::advectionDiffusion(dim, std::move(kappa), std::move(velocity)));
        p.mass1 = Number{1.0};
        return p;
    }
};

}  // namespace pde

#endif  // PDE_CORE_OPERATOR_HPP

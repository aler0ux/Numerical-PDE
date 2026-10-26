// Tests for the unified PDE representation (src/core).
//
// Build and run standalone until Phase 0 wires up CMake:
//   g++ -std=c++17 -Wall -Wextra -o test_operator tests/test_operator.cpp && ./test_operator

#include <cmath>
#include <iostream>
#include <limits>
#include <vector>

#include "../src/core/Operator.hpp"

using namespace pde;

namespace {

int failures = 0;

#define CHECK(cond)                                                        \
    do {                                                                   \
        if (!(cond)) {                                                     \
            std::cout << "FAIL line " << __LINE__ << ": " #cond "\n";      \
            ++failures;                                                    \
        }                                                                  \
    } while (false)

#define CHECK_NEAR(got, want)                                                          \
    do {                                                                               \
        const double g_ = (got), w_ = (want);                                          \
        if (!(std::fabs(g_ - w_) <= 1e-12)) {                                          \
            std::cout << "FAIL line " << __LINE__ << ": " #got " == " #want            \
                      << " (got " << g_ << ", want " << w_ << ")\n";                   \
            ++failures;                                                                 \
        }                                                                              \
    } while (false)

#define CHECK_EQ_STR(got, want)                                                        \
    do {                                                                               \
        const std::string g_ = (got), w_ = (want);                                     \
        if (g_ != w_) {                                                                \
            std::cout << "FAIL line " << __LINE__ << ": " #got                         \
                      << " (got \"" << g_ << "\", want \"" << w_ << "\")\n";           \
            ++failures;                                                                 \
        }                                                                              \
    } while (false)

#define CHECK_THROWS(expr)                                                 \
    do {                                                                   \
        bool threw_ = false;                                               \
        try {                                                              \
            (void)(expr);                                                  \
        } catch (const std::exception&) {                                  \
            threw_ = true;                                                 \
        }                                                                  \
        if (!threw_) {                                                     \
            std::cout << "FAIL line " << __LINE__                          \
                      << ": expected throw from " #expr "\n";              \
            ++failures;                                                    \
        }                                                                  \
    } while (false)

// Coefficient of d^alpha u, or NaN when the operator has no such term.
double coeffAt(const LinearOperator& L, const MultiIndex& alpha, const Env& env) {
    const Expr* c = L.coefficient(alpha);
    return c ? eval(*c, env) : std::numeric_limits<double>::quiet_NaN();
}

// -------------------------------------------------------------------------
void testExprEval() {
    double p[2] = {3.0, 4.0};
    const Env env{p, 2, 0.5};

    // The example from _scratch/test2.cpp: (x + 2) * 3 at x = 5.
    double q[1] = {5.0};
    CHECK_NEAR(eval((x() + 2.0) * 3.0, Env{q, 1, 0.0}), 21.0);

    CHECK_NEAR(eval(2.0 * x() + y() * y() - t(), env), 2 * 3.0 + 16.0 - 0.5);
    CHECK_NEAR(eval(pde::exp(-x()) * pde::sin(y()), env), std::exp(-3.0) * std::sin(4.0));
    CHECK_NEAR(eval(pde::pow(x(), 3.0), env), 27.0);
    CHECK_NEAR(eval(num(7.0)), 7.0);

    // Reaching for an axis the Env does not have is a bug, not a silent zero.
    CHECK_THROWS(eval(z(), env));
}

void testExprValueSemantics() {
    // The whole point of moving off unique_ptr: expressions are copyable, so a
    // PDE can be stored, passed by value, and handed to the GUI.
    const Expr e = 2.0 * x() + 1.0;
    std::vector<Expr> v;
    v.push_back(e);
    v.push_back(e);
    const std::vector<Expr> copy = v;

    double p[1] = {4.0};
    const Env env{p, 1, 0.0};
    CHECK_NEAR(eval(copy[0], env), 9.0);
    CHECK_NEAR(eval(copy[1], env), 9.0);
    CHECK_NEAR(eval(e, env), 9.0);
}

void testExprFolding() {
    CHECK(isZero(0.0 * x()));
    CHECK(isZero(x() * 0.0));
    CHECK(isZero(x() - x() * 1.0) == false);   // no algebraic simplification, only folding
    CHECK(std::holds_alternative<Coord>(1.0 * x()));
    CHECK(std::holds_alternative<Coord>(x() * 1.0));
    CHECK(std::holds_alternative<Coord>(x() + 0.0));
    CHECK(std::holds_alternative<Coord>(-(-x())));
    CHECK(isConstant(num(2.0) + num(3.0)));
    CHECK_NEAR(eval(num(2.0) + num(3.0)), 5.0);
    CHECK(isOne(pde::pow(x(), 0.0)));
}

// -------------------------------------------------------------------------
void testMultiIndex() {
    CHECK(MultiIndex::u().order() == 0);
    CHECK(MultiIndex::u().isIdentity());
    CHECK(MultiIndex::d(0).order() == 1);
    CHECK(MultiIndex::d2(0).order() == 2);
    CHECK(MultiIndex::d2(0, 1).order() == 2);

    CHECK_EQ_STR(MultiIndex::u().toString(), "u");
    CHECK_EQ_STR(MultiIndex::d(0).toString(), "du/dx");
    CHECK_EQ_STR(MultiIndex::d(1).toString(), "du/dy");
    CHECK_EQ_STR(MultiIndex::d2(0).toString(), "d2u/dx2");
    CHECK_EQ_STR(MultiIndex::d2(1).toString(), "d2u/dy2");
    CHECK_EQ_STR(MultiIndex::d2(0, 1).toString(), "d2u/dxdy");

    // Mixed partials canonicalise themselves -- the thing a (string, string)
    // factor pair could never do.
    CHECK(MultiIndex::d2(0, 1) == MultiIndex::d2(1, 0));
    CHECK(MultiIndex::d2(0) != MultiIndex::d2(0, 1));
    CHECK(MultiIndex::d(0) != MultiIndex::d(1));

    CHECK_THROWS(MultiIndex::d(kMaxDim));
    CHECK_THROWS(MultiIndex::d(-1));
}

// -------------------------------------------------------------------------
void testLinearOperator() {
    const Env env{nullptr, 0, 0.0};

    const LinearOperator lap = LinearOperator::laplacian(2);
    CHECK(lap.terms().size() == 2);
    CHECK(lap.maxOrder() == 2);
    CHECK_NEAR(coeffAt(lap, MultiIndex::d2(0), env), 1.0);
    CHECK_NEAR(coeffAt(lap, MultiIndex::d2(1), env), 1.0);
    CHECK(lap.coefficient(MultiIndex::d2(0, 1)) == nullptr);
    CHECK(lap.coefficient(MultiIndex::u()) == nullptr);

    // Same code path in 3D -- no new types, just a longer multi-index.
    const LinearOperator lap3 = LinearOperator::laplacian(3, num(0.5));
    CHECK(lap3.terms().size() == 3);
    CHECK_NEAR(coeffAt(lap3, MultiIndex::d2(2), env), 0.5);

    // Repeated multi-indices accumulate instead of duplicating.
    LinearOperator acc(2);
    acc.add(MultiIndex::d2(0), num(1.0));
    acc.add(MultiIndex::d2(0), num(2.0));
    CHECK(acc.terms().size() == 1);
    CHECK_NEAR(coeffAt(acc, MultiIndex::d2(0), env), 3.0);

    // ... and cancel out entirely when they sum to zero.
    acc.add(MultiIndex::d2(0), num(-3.0));
    CHECK(acc.terms().empty());

    // A zero coefficient never becomes a term (assembly should not walk it).
    LinearOperator z(2);
    z.add(MultiIndex::d(0), num(0.0));
    CHECK(z.terms().empty());

    // Terms are stored principal part first.
    LinearOperator mixed(2);
    mixed.add(MultiIndex::u(), num(1.0));
    mixed.add(MultiIndex::d(0), num(1.0));
    mixed.add(MultiIndex::d2(1), num(1.0));
    CHECK(mixed.terms().size() == 3);
    CHECK(mixed.terms()[0].alpha.order() == 2);
    CHECK(mixed.terms()[1].alpha.order() == 1);
    CHECK(mixed.terms()[2].alpha.isIdentity());

    // Using an axis the operator does not have is caught.
    LinearOperator oned(1);
    CHECK_THROWS(oned.add(MultiIndex::d(1), num(1.0)));
    CHECK_THROWS(LinearOperator(0));
    CHECK_THROWS(LinearOperator(kMaxDim + 1));

    // Variable coefficients.
    double p[2] = {2.0, 3.0};
    const Env at{p, 2, 0.0};
    LinearOperator vary(2);
    vary.add(MultiIndex::d2(0), 1.0 + x() * y());
    CHECK_NEAR(coeffAt(vary, MultiIndex::d2(0), at), 7.0);
}

// -------------------------------------------------------------------------
void testPDEBuilders() {
    const Env env{nullptr, 0, 0.0};

    const PDE poisson = PDE::poisson(2, num(3.0));
    CHECK(poisson.timeKind() == TimeKind::Steady);
    CHECK(poisson.isSteady());
    // u_xx + u_yy = 3  is stored as  0 = u_xx + u_yy - 3.
    CHECK_NEAR(eval(poisson.spatial.source(), env), -3.0);

    const PDE heat = PDE::heat(2, num(0.1));
    CHECK(heat.timeKind() == TimeKind::FirstOrder);
    CHECK_NEAR(eval(heat.mass1, env), 1.0);
    CHECK(isZero(heat.mass2));
    CHECK(heat.spatial.terms().size() == 2);
    CHECK_NEAR(coeffAt(heat.spatial, MultiIndex::d2(0), env), 0.1);

    const PDE wave = PDE::wave(2, num(2.0));
    CHECK(wave.timeKind() == TimeKind::SecondOrder);
    CHECK_NEAR(eval(wave.mass2, env), 1.0);
    CHECK_NEAR(coeffAt(wave.spatial, MultiIndex::d2(0), env), 4.0);  // speed^2, folded

    const PDE adv = PDE::advectionDiffusion(2, num(0.01), {num(1.0), num(-2.0)});
    CHECK(adv.timeKind() == TimeKind::FirstOrder);
    CHECK_NEAR(coeffAt(adv.spatial, MultiIndex::d2(1), env), 0.01);
    CHECK_NEAR(coeffAt(adv.spatial, MultiIndex::d(0), env), -1.0);
    CHECK_NEAR(coeffAt(adv.spatial, MultiIndex::d(1), env), 2.0);
    CHECK_THROWS(PDE::advectionDiffusion(2, num(1.0), {num(1.0)}));

    const PDE helm = PDE::steady(LinearOperator::helmholtz(2, num(9.0)));
    CHECK(helm.isSteady());
    CHECK_NEAR(coeffAt(helm.spatial, MultiIndex::u(), env), 9.0);
}

// -------------------------------------------------------------------------
// The two structs this representation replaces, expressed in the new form.
void testSubsumesOldStructs() {
    const Env env{nullptr, 0, 0.0};

    // GeneralizedSecondOrder was:
    //   A u_xx + B u_xy + C u_yy + D u_x + E u_y + F u = R(x, y)
    // Six named fields become six terms, and the source flips sign because the
    // canonical form keeps everything on the right.
    const double A = 1.0, B = 0.5, C = 2.0, D = 3.0, E = 4.0, F = 5.0, R = 6.0;
    LinearOperator second(2);
    second.add(MultiIndex::d2(0), num(A));
    second.add(MultiIndex::d2(0, 1), num(B));
    second.add(MultiIndex::d2(1), num(C));
    second.add(MultiIndex::d(0), num(D));
    second.add(MultiIndex::d(1), num(E));
    second.add(MultiIndex::u(), num(F));
    second.setSource(num(-R));

    const PDE p2 = PDE::steady(second);
    CHECK(p2.spatial.terms().size() == 6);
    CHECK(p2.isSteady());
    CHECK_NEAR(coeffAt(p2.spatial, MultiIndex::d2(0, 1), env), B);
    CHECK_NEAR(eval(p2.spatial.source(), env), -R);

    // GeneralizedFirstOrder was:
    //   a u_x + b u_y + c u_t + d u = r(x, y, t)
    // which rearranges to   c u_t = -a u_x - b u_y - d u + r.
    // The old struct could not express this alongside second-order terms; here
    // it is the same type, and second-order terms could simply be added.
    const double a = 2.0, b = 3.0, c = 4.0, d = 5.0, r = 6.0;
    LinearOperator firstSpatial(2);
    firstSpatial.add(MultiIndex::d(0), num(-a));
    firstSpatial.add(MultiIndex::d(1), num(-b));
    firstSpatial.add(MultiIndex::u(), num(-d));
    firstSpatial.setSource(num(r));

    PDE p1 = PDE::steady(firstSpatial);
    p1.mass1 = num(c);
    CHECK(p1.timeKind() == TimeKind::FirstOrder);
    CHECK(p1.spatial.maxOrder() == 1);
    CHECK_NEAR(eval(p1.mass1, env), c);
    CHECK_NEAR(coeffAt(p1.spatial, MultiIndex::d(0), env), -a);
    CHECK_NEAR(eval(p1.spatial.source(), env), r);

    // And the thing neither old struct could do: time derivative plus a full
    // second-order anisotropic spatial part, with variable coefficients.
    PDE hybrid = PDE::steady(LinearOperator::laplacian(2, 1.0 + pde::sin(x())));
    hybrid.spatial.add(MultiIndex::d2(0, 1), 0.5 * y());
    hybrid.spatial.add(MultiIndex::d(0), -pde::exp(-t()));
    hybrid.mass1 = num(1.0);
    CHECK(hybrid.timeKind() == TimeKind::FirstOrder);
    CHECK(hybrid.spatial.terms().size() == 4);
}

// -------------------------------------------------------------------------
// Printing is a debugging tool for assembly, so its output is pinned. The
// bracketing cases below are all ones that read as a different equation if
// they come out wrong.
void testPrinting() {
    // A sum used as a coefficient must bracket, or it swallows the derivative:
    // `1 + 0.5*sin(x)*d2u/dx2` is a different expression entirely.
    LinearOperator sum(2);
    sum.add(MultiIndex::d2(0), 1.0 + 0.5 * pde::sin(x()));
    CHECK_EQ_STR(sum.toString(), "(1 + 0.5*sin(x))*d2u/dx2");

    // A product coefficient needs no brackets.
    LinearOperator prod(2);
    prod.add(MultiIndex::d2(0, 1), 0.25 * y());
    CHECK_EQ_STR(prod.toString(), "0.25*y*d2u/dxdy");

    // Negative coefficients and sources print as a subtraction.
    CHECK_EQ_STR(PDE::poisson(2, num(3.0)).toString(), "0 = d2u/dx2 + d2u/dy2 - 3");
    CHECK_EQ_STR(PDE::heat(2, num(0.1)).toString(), "du/dt = 0.1*d2u/dx2 + 0.1*d2u/dy2");
    CHECK_EQ_STR(PDE::wave(2, num(2.0)).toString(), "d2u/dt2 = 4*d2u/dx2 + 4*d2u/dy2");
    CHECK_EQ_STR(PDE::advectionDiffusion(1, num(0.01), {num(1.0)}).toString(),
                 "du/dt = 0.01*d2u/dx2 - du/dx");

    // A leading negative term keeps its sign without a stray separator.
    LinearOperator lead(1);
    lead.add(MultiIndex::d2(0), num(-1.0));
    CHECK_EQ_STR(lead.toString(), "-d2u/dx2");

    // A negated sum as a source must keep its brackets.
    LinearOperator negSource(2);
    negSource.add(MultiIndex::d2(0), num(1.0));
    negSource.setSource(-(x() + y()));
    CHECK_EQ_STR(negSource.toString(), "d2u/dx2 - (x + y)");

    // Textbook term order, mixed partial included.
    LinearOperator full(2);
    full.add(MultiIndex::d(1), num(1.0));
    full.add(MultiIndex::u(), num(1.0));
    full.add(MultiIndex::d2(1), num(1.0));
    full.add(MultiIndex::d2(0, 1), num(1.0));
    full.add(MultiIndex::d(0), num(1.0));
    full.add(MultiIndex::d2(0), num(1.0));
    CHECK_EQ_STR(full.toString(), "d2u/dx2 + d2u/dxdy + d2u/dy2 + du/dx + du/dy + u");

    LinearOperator empty(1);
    CHECK_EQ_STR(empty.toString(), "0");
}

void showOperators() {
    std::cout << "\nCanonical form, printed:\n";
    std::cout << "  poisson(2)      " << PDE::poisson(2, num(3.0)).toString() << "\n";
    std::cout << "  heat(2)         " << PDE::heat(2, num(0.1)).toString() << "\n";
    std::cout << "  wave(2)         " << PDE::wave(2, num(2.0)).toString() << "\n";
    std::cout << "  wave(3)         " << PDE::wave(3, num(1.0)).toString() << "\n";
    std::cout << "  adv-diff(2)     "
              << PDE::advectionDiffusion(2, num(0.01), {num(1.0), num(-2.0)}).toString() << "\n";
    std::cout << "  helmholtz(2)    "
              << PDE::steady(LinearOperator::helmholtz(2, num(9.0))).toString() << "\n";

    PDE varying = PDE::heat(2, 1.0 + 0.5 * pde::sin(x()));
    varying.spatial.add(MultiIndex::d2(0, 1), 0.25 * y());
    std::cout << "  anisotropic     " << varying.toString() << "\n";
}

}  // namespace

int main() {
    testExprEval();
    testExprValueSemantics();
    testExprFolding();
    testMultiIndex();
    testLinearOperator();
    testPDEBuilders();
    testSubsumesOldStructs();
    testPrinting();

    if (failures == 0) {
        std::cout << "all tests passed\n";
        showOperators();
        return 0;
    }
    std::cout << failures << " check(s) failed\n";
    return 1;
}

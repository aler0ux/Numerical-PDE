#ifndef PDE_CORE_EXPR_HPP
#define PDE_CORE_EXPR_HPP

// Scalar expressions over space and time.
//
// An Expr is what every user-supplied scalar in the library is made of: PDE
// coefficients, source terms, boundary data, initial conditions. It is a small
// immutable tree with value semantics -- copying an Expr is a pointer bump, and
// an Expr can be stored in containers, which is what lets a whole PDE be a
// copyable value.

#include <cmath>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <utility>
#include <variant>

namespace pde {

// ---------------------------------------------------------------------------
// Evaluation context: one point in space-time.
//
// Coordinates are addressed by index rather than by name so that eval() costs
// an array load instead of a string comparison -- it runs once per matrix entry
// during assembly.
// ---------------------------------------------------------------------------
struct Env {
    const double* x = nullptr;  // spatial coordinates, x[0 .. dim-1]
    int           dim = 0;
    double        t = 0.0;
};

// ---------------------------------------------------------------------------
// Nodes
// ---------------------------------------------------------------------------
struct Binary;
struct Unary;

struct Number { double val = 0.0; };  // a constant
struct Coord  { int axis = 0; };      // the spatial coordinate x_axis
struct Time   {};                     // t

using Expr = std::variant<Number, Coord, Time,
                          std::shared_ptr<Binary>,
                          std::shared_ptr<Unary>>;

enum class BinOp { Add, Sub, Mul, Div, Pow };
enum class UnFn { Neg, Sin, Cos, Exp, Log, Sqrt, Tanh, Abs, Floor };

struct Binary {
    BinOp op;
    Expr  lhs, rhs;
};

struct Unary {
    UnFn fn;
    Expr arg;
};

// ---------------------------------------------------------------------------
// Scalar kernels, shared by eval() and by constant folding
// ---------------------------------------------------------------------------
inline double applyBinary(BinOp op, double a, double b) {
    switch (op) {
        case BinOp::Add: return a + b;
        case BinOp::Sub: return a - b;
        case BinOp::Mul: return a * b;
        case BinOp::Div: return a / b;
        case BinOp::Pow: return std::pow(a, b);
    }
    throw std::logic_error("pde::applyBinary: unhandled BinOp");
}

inline double applyUnary(UnFn fn, double a) {
    switch (fn) {
        case UnFn::Neg:  return -a;
        case UnFn::Sin:  return std::sin(a);
        case UnFn::Cos:  return std::cos(a);
        case UnFn::Exp:  return std::exp(a);
        case UnFn::Log:  return std::log(a);
        case UnFn::Sqrt: return std::sqrt(a);
        case UnFn::Tanh: return std::tanh(a);
        case UnFn::Abs:  return std::fabs(a);
        case UnFn::Floor: return std::floor(a);
    }
    throw std::logic_error("pde::applyUnary: unhandled UnFn");
}

// ---------------------------------------------------------------------------
// Constant inspection
// ---------------------------------------------------------------------------
inline bool asNumber(const Expr& e, double& out) {
    if (const Number* n = std::get_if<Number>(&e)) {
        out = n->val;
        return true;
    }
    return false;
}

inline bool isConstant(const Expr& e) { return std::holds_alternative<Number>(e); }

inline bool isZero(const Expr& e) {
    double v = 0.0;
    return asNumber(e, v) && v == 0.0;
}

inline bool isOne(const Expr& e) {
    double v = 0.0;
    return asNumber(e, v) && v == 1.0;
}

// ---------------------------------------------------------------------------
// Construction, with enough folding that generated operators stay readable
// (a Laplacian built from a unit coefficient should not carry `1 *` nodes).
// ---------------------------------------------------------------------------
inline Expr makeUnary(UnFn fn, Expr a) {
    double v = 0.0;
    if (asNumber(a, v)) return Number{applyUnary(fn, v)};
    if (fn == UnFn::Neg) {
        if (auto* u = std::get_if<std::shared_ptr<Unary>>(&a)) {
            if ((*u)->fn == UnFn::Neg) return (*u)->arg;  // -(-e) == e
        }
    }
    return std::make_shared<Unary>(Unary{fn, std::move(a)});
}

inline Expr makeBinary(BinOp op, Expr a, Expr b) {
    double va = 0.0, vb = 0.0;
    const bool ca = asNumber(a, va);
    const bool cb = asNumber(b, vb);
    if (ca && cb) return Number{applyBinary(op, va, vb)};

    switch (op) {
        case BinOp::Add:
            if (isZero(a)) return b;
            if (isZero(b)) return a;
            break;
        case BinOp::Sub:
            if (isZero(b)) return a;
            if (isZero(a)) return makeUnary(UnFn::Neg, std::move(b));
            break;
        case BinOp::Mul:
            if (isZero(a) || isZero(b)) return Number{0.0};
            if (isOne(a)) return b;
            if (isOne(b)) return a;
            break;
        case BinOp::Div:
            if (isZero(a)) return Number{0.0};
            if (isOne(b)) return a;
            break;
        case BinOp::Pow:
            if (isZero(b)) return Number{1.0};
            if (isOne(b)) return a;
            break;
    }
    return std::make_shared<Binary>(Binary{op, std::move(a), std::move(b)});
}

// ---------------------------------------------------------------------------
// Leaves
// ---------------------------------------------------------------------------
inline Expr num(double v) { return Number{v}; }
inline Expr axis(int i) { return Coord{i}; }
inline Expr x() { return Coord{0}; }
inline Expr y() { return Coord{1}; }
inline Expr z() { return Coord{2}; }
inline Expr t() { return Time{}; }

// ---------------------------------------------------------------------------
// Operators. The double-mixed overloads exist because double -> Number -> Expr
// is two user-defined conversions and so never happens implicitly.
// ---------------------------------------------------------------------------
inline Expr operator+(Expr a, Expr b) { return makeBinary(BinOp::Add, std::move(a), std::move(b)); }
inline Expr operator-(Expr a, Expr b) { return makeBinary(BinOp::Sub, std::move(a), std::move(b)); }
inline Expr operator*(Expr a, Expr b) { return makeBinary(BinOp::Mul, std::move(a), std::move(b)); }
inline Expr operator/(Expr a, Expr b) { return makeBinary(BinOp::Div, std::move(a), std::move(b)); }
inline Expr operator-(Expr a) { return makeUnary(UnFn::Neg, std::move(a)); }

inline Expr operator+(Expr a, double b) { return std::move(a) + num(b); }
inline Expr operator+(double a, Expr b) { return num(a) + std::move(b); }
inline Expr operator-(Expr a, double b) { return std::move(a) - num(b); }
inline Expr operator-(double a, Expr b) { return num(a) - std::move(b); }
inline Expr operator*(Expr a, double b) { return std::move(a) * num(b); }
inline Expr operator*(double a, Expr b) { return num(a) * std::move(b); }
inline Expr operator/(Expr a, double b) { return std::move(a) / num(b); }
inline Expr operator/(double a, Expr b) { return num(a) / std::move(b); }

inline Expr pow(Expr a, Expr b) { return makeBinary(BinOp::Pow, std::move(a), std::move(b)); }
inline Expr pow(Expr a, double b) { return pow(std::move(a), num(b)); }
inline Expr sin(Expr a) { return makeUnary(UnFn::Sin, std::move(a)); }
inline Expr cos(Expr a) { return makeUnary(UnFn::Cos, std::move(a)); }
inline Expr exp(Expr a) { return makeUnary(UnFn::Exp, std::move(a)); }
inline Expr log(Expr a) { return makeUnary(UnFn::Log, std::move(a)); }
inline Expr sqrt(Expr a) { return makeUnary(UnFn::Sqrt, std::move(a)); }
inline Expr tanh(Expr a) { return makeUnary(UnFn::Tanh, std::move(a)); }
inline Expr abs(Expr a) { return makeUnary(UnFn::Abs, std::move(a)); }
inline Expr floor(Expr a) { return makeUnary(UnFn::Floor, std::move(a)); }

// ---------------------------------------------------------------------------
// Evaluation
// ---------------------------------------------------------------------------
inline double eval(const Expr& e, const Env& env) {
    return std::visit(
        [&env](const auto& node) -> double {
            using T = std::decay_t<decltype(node)>;
            if constexpr (std::is_same_v<T, Number>) {
                return node.val;
            } else if constexpr (std::is_same_v<T, Coord>) {
                if (!env.x || node.axis < 0 || node.axis >= env.dim)
                    throw std::out_of_range("pde::eval: coordinate axis outside Env");
                return env.x[node.axis];
            } else if constexpr (std::is_same_v<T, Time>) {
                return env.t;
            } else if constexpr (std::is_same_v<T, std::shared_ptr<Binary>>) {
                return applyBinary(node->op, eval(node->lhs, env), eval(node->rhs, env));
            } else {
                return applyUnary(node->fn, eval(node->arg, env));
            }
        },
        e);
}

// Convenience for constant expressions and for 1D/2D call sites.
inline double eval(const Expr& e, double t_value = 0.0) {
    return eval(e, Env{nullptr, 0, t_value});
}

// ---------------------------------------------------------------------------
// Printing (debugging assembly without being able to print the operator is
// miserable, so this exists from day one)
// ---------------------------------------------------------------------------
inline const char* axisName(int i) {
    static const char* names[] = {"x", "y", "z", "w"};
    return (i >= 0 && i < 4) ? names[i] : "?";
}

inline std::string formatNumber(double v) {
    std::ostringstream os;
    os << v;
    return os.str();
}

// Binding strength, so printing only parenthesises where it must.
inline int precedence(BinOp op) {
    switch (op) {
        case BinOp::Add:
        case BinOp::Sub: return 1;
        case BinOp::Mul:
        case BinOp::Div: return 2;
        case BinOp::Pow: return 3;
    }
    return 0;
}

// `parentPrec` is the binding strength of the context this expression sits in.
inline std::string toStringAt(const Expr& e, int parentPrec) {
    return std::visit(
        [parentPrec](const auto& node) -> std::string {
            using T = std::decay_t<decltype(node)>;
            if constexpr (std::is_same_v<T, Number>) {
                return formatNumber(node.val);
            } else if constexpr (std::is_same_v<T, Coord>) {
                return axisName(node.axis);
            } else if constexpr (std::is_same_v<T, Time>) {
                return "t";
            } else if constexpr (std::is_same_v<T, std::shared_ptr<Binary>>) {
                const char* sym = "?";
                // Subtraction and division do not reassociate on the right, so
                // a right child of equal strength still needs brackets.
                int rightBias = 0;
                switch (node->op) {
                    case BinOp::Add: sym = " + "; break;
                    case BinOp::Sub: sym = " - "; rightBias = 1; break;
                    case BinOp::Mul: sym = "*"; break;
                    case BinOp::Div: sym = "/"; rightBias = 1; break;
                    case BinOp::Pow: sym = "^"; break;
                }
                const int prec = precedence(node->op);
                const std::string s = toStringAt(node->lhs, prec) + sym +
                                      toStringAt(node->rhs, prec + rightBias);
                return prec < parentPrec ? "(" + s + ")" : s;
            } else {
                const char* name = "?";
                switch (node->fn) {
                    case UnFn::Neg: {
                        // Treated as loosely binding so that it brackets inside
                        // a product: 2*(-x) rather than 2*-x.
                        const std::string s = "-" + toStringAt(node->arg, 2);
                        return 1 < parentPrec ? "(" + s + ")" : s;
                    }
                    case UnFn::Sin:  name = "sin"; break;
                    case UnFn::Cos:  name = "cos"; break;
                    case UnFn::Exp:  name = "exp"; break;
                    case UnFn::Log:  name = "log"; break;
                    case UnFn::Sqrt: name = "sqrt"; break;
                    case UnFn::Tanh: name = "tanh"; break;
                    case UnFn::Abs:  name = "abs"; break;
                    case UnFn::Floor: name = "floor"; break;
                }
                return name + ("(" + toStringAt(node->arg, 0) + ")");
            }
        },
        e);
}

inline std::string toString(const Expr& e) { return toStringAt(e, 0); }

// Peel a leading minus off an expression so that sums can be printed as
// `a - b` instead of `a + -b`. Returns true when a sign was removed.
inline bool splitSign(const Expr& e, Expr& magnitude) {
    double v = 0.0;
    if (asNumber(e, v) && v < 0.0) {
        magnitude = Number{-v};
        return true;
    }
    if (auto* u = std::get_if<std::shared_ptr<Unary>>(&e)) {
        if ((*u)->fn == UnFn::Neg) {
            magnitude = (*u)->arg;
            return true;
        }
    }
    magnitude = e;
    return false;
}

}  // namespace pde

#endif  // PDE_CORE_EXPR_HPP

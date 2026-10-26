#ifndef PDE_COMPARE_SOURCE_FUNCTIONS_HPP
#define PDE_COMPARE_SOURCE_FUNCTIONS_HPP

// A small, selectable library of right-hand sides f(x,y) for
//
//     -Laplacian(u) = f    on the domain,     u = 0 on the boundary
//
// Each one is a `pde::Expr` -- the same symbolic expression type the rest
// of the library (core/Expr.hpp, core/Operator.hpp) uses for coefficients
// and sources, so picking a different f here is exactly the "configurable
// f" the project's PDE framework is already built to support, not a
// separate ad hoc mechanism.

#include <algorithm>
#include <stdexcept>
#include <string>
#include <vector>

#include "../core/Expr.hpp"

namespace cmp {

// Local shorthand for pi, used when building sin(pi*x)-style source terms
// below. Kept as our own named constant (rather than <cmath>'s M_PI, which
// isn't standard C++ and isn't guaranteed to be defined on every compiler,
// e.g. MSVC without _USE_MATH_DEFINES) so it is portable across the
// Windows/Linux builds this project targets.
inline constexpr double kPi = 3.14159265358979323846;

struct SourceFunction {
    std::string name;
    std::string description;
    pde::Expr expr;  // function of x = Coord{0}, y = Coord{1}
};

inline const std::vector<SourceFunction>& sourceFunctionLibrary() {
    using namespace pde;
    static const std::vector<SourceFunction> lib = {
        {"constant", "f = 1 (uniform load)", num(1.0)},
        {"negative constant", "f = -9.8 (gravity)", num(-9.8)},
        {"sheer", "f = floor(x+0.5) + floor(y+0.5)", floor(num(0.5) + x()) + floor(num(0.5) + y())},
        {"sine", "f = sin(pi x) sin(pi y)", sin(num(kPi) * x()) * sin(num(kPi) * y())},
        {"bump", "f = exp(-((x-0.5)^2+(y-0.5)^2)/0.02) (a Gaussian bump at the centre)",
         exp(-((x() - 0.5) * (x() - 0.5) + (y() - 0.5) * (y() - 0.5)) / 0.02)},
        {"linear", "f = x + y", x() + y()},
        {"checker", "f = sin(2 pi x) * sin(3 pi y)",
         sin(num(2.0 * kPi) * x()) * sin(num(3.0 * kPi) * y())},
    };
    return lib;
}

inline const SourceFunction& sourceFunctionFromName(const std::string& nameIn) {
    std::string name = nameIn;
    std::transform(name.begin(), name.end(), name.begin(), [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
    for (const auto& f : sourceFunctionLibrary())
        if (f.name == name) return f;
    throw std::invalid_argument("unknown source function '" + nameIn + "'");
}

inline double evalF(const pde::Expr& f, double x, double y) {
    const double coords[2] = {x, y};
    return pde::eval(f, pde::Env{coords, 2, 0.0});
}

}  // namespace cmp

#endif  // PDE_COMPARE_SOURCE_FUNCTIONS_HPP

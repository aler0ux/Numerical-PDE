#ifndef PDE_DOMAIN2D_HPP
#define PDE_DOMAIN2D_HPP

// 2D CSG domain: boolean add (union) / subtract (difference) of simple
// shapes (axis-aligned rectangles, circles), evaluated as a left-to-right
// fold over the recorded operations -- not a general CSG tree, so order
// matters: a subtract only removes area already covered by an earlier add.
//
// Previously split across Domain2d.hpp (declarations) and Domain2d.cpp
// (out-of-line definitions); merged into one header-only file here, to
// match core/Expr.hpp and core/Operator.hpp's style and because nothing in
// the project currently compiles Domain2d.cpp as its own translation unit
// (Domain2d.cpp is now just a pointer to this file). This also fixes the
// bug inventory recorded in docs/ROADMAP.md:
//   - inDomain() computed `is_inside` but fell off the end without
//     returning it (undefined behaviour).
//   - Circle::contains compared squared distance against `radius` instead
//     of `radius*radius`.
//   - Circle::on_boundary used exact `==` on doubles, effectively never
//     true; now uses the same epsilon tolerance Rectangle::on_boundary
//     already used.
//   - The subtract case in inDomain() computed `is_inside && inOperation`
//     (an intersection) instead of `is_inside && !inOperation`.
//   - onBoundary() computed `onOperationBoundary` and never read it, so
//     every branch only ever assigned `false` or left the value
//     unchanged -- it always returned false. See that method's own
//     comment for how it is fixed here.
//   - Both Domain1d.hpp and (the old) Domain2d.hpp used the include guard
//     `DOMAIN_HPP` and both defined `pde::EPSILON`, so including both in
//     one translation unit silently dropped whichever was included
//     second. This file now uses its own guard and its own
//     epsilon constant, `kDomain2dEpsilon`, so it no longer depends on
//     (or collides with) anything Domain1d.hpp defines.

#include <algorithm>
#include <cmath>
#include <variant>
#include <vector>

namespace pde {

inline constexpr double kDomain2dEpsilon = 1e-9;

struct Rectangle {
    double x_min, x_max;
    double y_min, y_max;

    bool contains(double x, double y) const {
        return x >= x_min && x <= x_max && y >= y_min && y <= y_max;
    }
    bool on_boundary(double x, double y) const {
        return contains(x, y) &&
               (std::abs(x - x_min) < kDomain2dEpsilon || std::abs(x - x_max) < kDomain2dEpsilon ||
                std::abs(y - y_min) < kDomain2dEpsilon || std::abs(y - y_max) < kDomain2dEpsilon);
    }
};

struct Circle {
    double x_center, y_center;
    double radius;

    bool contains(double x, double y) const {
        const double dx = x - x_center;
        const double dy = y - y_center;
        return dx * dx + dy * dy <= radius * radius;
    }
    bool on_boundary(double x, double y) const {
        const double dx = x - x_center;
        const double dy = y - y_center;
        return std::abs(std::sqrt(dx * dx + dy * dy) - radius) < kDomain2dEpsilon;
    }
};

// TODO: add more shapes -- must be included in Shape below too.
using Shape = std::variant<Rectangle, Circle>;

enum class OpType { Add, Subtract };

struct DomainOp {
    OpType type;
    Shape shape;
};

class Domain2d {
public:
    void addShape(Shape s) { domainOperations.push_back({OpType::Add, s}); }
    void removeShape(Shape s) { domainOperations.push_back({OpType::Subtract, s}); }

    // Left-to-right fold of the recorded operations: an Add unions the
    // shape's interior in, a Subtract removes it from whatever has been
    // accumulated so far. A Subtract before any Add therefore does
    // nothing (there is nothing yet to remove from) -- shapes should be
    // added in the order a CSG expression would be read left to right.
    bool inDomain(double x, double y) const {
        bool is_inside = false;
        for (const auto& op : domainOperations) {
            const bool inOperation =
                std::visit([x, y](const auto& shape) { return shape.contains(x, y); }, op.shape);
            if (op.type == OpType::Add) {
                is_inside = is_inside || inOperation;
            } else {
                is_inside = is_inside && !inOperation;
            }
        }
        return is_inside;
    }

    // A point is "on the boundary" of the CSG result if it is itself in
    // the (closed) domain and an arbitrarily small nudge in some
    // direction leaves it. This sidesteps tracking which primitive's
    // surface "survives" the remaining operations -- the general, harder
    // version of this function, and the one piece of the geometry layer
    // the roadmap calls out as needing a signed-distance formulation to
    // get right -- at the cost of needing a probe distance much smaller
    // than whatever grid spacing this is tested against. That is exactly
    // how this class is used in practice (src/compare/Grid2d.hpp's own
    // node classification calls only inDomain(), for precisely this
    // reason), so the simplification is appropriate here even though it
    // is not the fully general answer.
    bool onBoundary(double x, double y) const {
        if (!inDomain(x, y)) return false;
        constexpr double kProbe = 1e-6;
        return !inDomain(x + kProbe, y) || !inDomain(x - kProbe, y) || !inDomain(x, y + kProbe) ||
               !inDomain(x, y - kProbe);
    }

private:
    std::vector<DomainOp> domainOperations;
};

}  // namespace pde

#endif  // PDE_DOMAIN2D_HPP

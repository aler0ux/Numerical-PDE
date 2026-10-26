#ifndef PDE_COMPARE_DOMAINS_HPP
#define PDE_COMPARE_DOMAINS_HPP

// Four preset 2D domains, each built as a small CSG expression -- unions
// and differences of rectangles and circles -- using the project's own
// `pde::Domain2d` (src/Domain2d.hpp), exactly as it was written: add a
// shape to union it in, remove a shape to subtract it from whatever has
// been added so far.
//
// An earlier version of this file hand-wrote four point-in-domain
// predicates instead, because Domain2d.hpp had several open correctness
// bugs at the time (see docs/ROADMAP.md's bug inventory: a missing return
// in inDomain(), a squared-vs-unsquared radius comparison in Circle, an
// exact `==` comparison on doubles, and a subtract that behaved as an
// intersect). Those are now fixed directly in Domain2d.hpp, so this file
// goes back to using it -- the whole point of having a CSG domain type is
// not needing a bespoke predicate per shape.

#include <algorithm>
#include <stdexcept>
#include <string>

#include "../Domain2d.hpp"

namespace cmp {

enum class DomainKind { Square, LShape, Plus, Circle };

inline std::string domainName(DomainKind k) {
    switch (k) {
        case DomainKind::Square: return "square";
        case DomainKind::LShape: return "lshape";
        case DomainKind::Plus: return "plus";
        case DomainKind::Circle: return "circle";
    }
    return "?";
}

inline DomainKind domainFromName(const std::string& nameIn) {
    std::string name = nameIn;
    std::transform(name.begin(), name.end(), name.begin(), [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
    if (name == "square" || name == "sq") return DomainKind::Square;
    if (name == "lshape" || name == "l" || name == "l-shape") return DomainKind::LShape;
    if (name == "plus" || name == "cross" || name == "+") return DomainKind::Plus;
    if (name == "circle" || name == "disk" || name == "c") return DomainKind::Circle;
    throw std::invalid_argument("unknown domain '" + nameIn + "' (expected square, lshape, plus, circle)");
}

// Build the CSG domain for one of the four presets. All four live inside
// [0,1] x [0,1].
inline pde::Domain2d buildDomain2d(DomainKind kind) {
    pde::Domain2d d;
    switch (kind) {
        case DomainKind::Square: {
            d.addShape(pde::Rectangle{0.0, 1.0, 0.0, 1.0});
            break;
        }
        case DomainKind::LShape: {
            // Unit square, add; top-right quadrant, subtract -- a single
            // add-then-remove, read left to right exactly like the CSG
            // expression it represents. The removed rectangle is nudged
            // in by a hair so the notch's own edge (x=0.5 or y=0.5)
            // stays *in* the domain rather than being shaved off by
            // Rectangle::contains' closed (>=, <=) interval test -- that
            // edge is meant to be part of the boundary, not excluded
            // from it.
            constexpr double kNudge = 1e-6;
            d.addShape(pde::Rectangle{0.0, 1.0, 0.0, 1.0});
            d.removeShape(pde::Rectangle{0.5 + kNudge, 1.0, 0.5 + kNudge, 1.0});
            break;
        }
        case DomainKind::Plus: {
            // Union of a horizontal and a vertical bar through the centre.
            d.addShape(pde::Rectangle{0.2, 0.8, 0.4, 0.6});
            d.addShape(pde::Rectangle{0.4, 0.6, 0.2, 0.8});
            break;
        }
        case DomainKind::Circle: {
            d.addShape(pde::Circle{0.5, 0.5, 0.45});
            break;
        }
    }
    return d;
}

}  // namespace cmp

#endif  // PDE_COMPARE_DOMAINS_HPP

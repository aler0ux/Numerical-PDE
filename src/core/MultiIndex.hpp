#ifndef PDE_CORE_MULTIINDEX_HPP
#define PDE_CORE_MULTIINDEX_HPP

// A derivative multi-index: which derivative of u a term applies to.
//
// alpha[i] is the order of differentiation along spatial axis i, so in 2D:
//
//     u            {0, 0}      order 0
//     du/dx        {1, 0}      order 1
//     du/dy        {0, 1}      order 1
//     d2u/dx2      {2, 0}      order 2
//     d2u/dxdy     {1, 1}      order 2
//     d2u/dy2      {0, 2}      order 2
//
// Keying terms this way instead of by name (A, B, C, ...) is what makes one
// code path cover every linear PDE: the representation does not care how many
// dimensions there are, mixed partials canonicalise themselves (d/dxdy and
// d/dydx are the same multi-index), and each multi-index maps directly onto a
// finite-difference stencil via a tensor product of 1D stencils.

#include <array>
#include <cstdint>
#include <stdexcept>
#include <string>

namespace pde {

// Largest number of spatial dimensions supported. Time is tracked separately
// (see PDE in Operator.hpp), so this counts x, y, z, ... only.
inline constexpr int kMaxDim = 4;

class MultiIndex {
public:
    MultiIndex() : a_{} {}

    // u itself -- the zeroth derivative.
    static MultiIndex u() { return MultiIndex{}; }

    // First derivative along one axis.
    static MultiIndex d(int ax) {
        MultiIndex m;
        m.bump(ax);
        return m;
    }

    // Second derivative. d2(0, 1) is the mixed partial; d2(0) is d2u/dx2.
    static MultiIndex d2(int ax1, int ax2) {
        MultiIndex m;
        m.bump(ax1);
        m.bump(ax2);
        return m;
    }
    static MultiIndex d2(int ax) { return d2(ax, ax); }

    int operator[](int ax) const {
        check(ax);
        return a_[static_cast<size_t>(ax)];
    }

    int order() const {
        int s = 0;
        for (uint8_t v : a_) s += v;
        return s;
    }

    bool isIdentity() const { return order() == 0; }

    friend bool operator==(const MultiIndex& p, const MultiIndex& q) { return p.a_ == q.a_; }
    friend bool operator!=(const MultiIndex& p, const MultiIndex& q) { return !(p == q); }

    // Lexicographic total order, so a MultiIndex can key a std::map or be sorted.
    friend bool operator<(const MultiIndex& p, const MultiIndex& q) { return p.a_ < q.a_; }

    // "u", "du/dx", "d2u/dx2", "d2u/dxdy", ...
    std::string toString() const {
        const int n = order();
        if (n == 0) return "u";
        std::string s = (n == 1) ? "du/" : ("d" + std::to_string(n) + "u/");
        for (int i = 0; i < kMaxDim; ++i) {
            const int k = a_[static_cast<size_t>(i)];
            if (k == 0) continue;
            s += "d";
            s += axisLabel(i);
            if (k > 1) s += std::to_string(k);
        }
        return s;
    }

private:
    void bump(int ax) {
        check(ax);
        ++a_[static_cast<size_t>(ax)];
    }

    static void check(int ax) {
        if (ax < 0 || ax >= kMaxDim)
            throw std::out_of_range("pde::MultiIndex: axis " + std::to_string(ax) +
                                    " outside [0, " + std::to_string(kMaxDim) + ")");
    }

    static const char* axisLabel(int i) {
        static const char* names[] = {"x", "y", "z", "w"};
        return (i >= 0 && i < 4) ? names[i] : "?";
    }

    std::array<uint8_t, kMaxDim> a_;
};

// Ordering used for display and canonical storage: highest-order terms first
// (the principal part), and within one order, earlier axes first. That yields
// the conventional textbook reading order -- u_xx, u_xy, u_yy, u_x, u_y, u --
// which needs descending lexicographic order, since a larger count in axis 0
// makes an index lexicographically greater.
struct PrincipalFirst {
    bool operator()(const MultiIndex& p, const MultiIndex& q) const {
        if (p.order() != q.order()) return p.order() > q.order();
        return q < p;
    }
};

}  // namespace pde

#endif  // PDE_CORE_MULTIINDEX_HPP

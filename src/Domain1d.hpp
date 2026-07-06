#ifndef DOMAIN_HPP
#define DOMAIN_HPP
#include <vector>
#include <variant>
#include <cmath>
namespace pde {
    constexpr double EPSILON = 1e-9;
    struct Interval {
        double x_min,x_max;
        bool contains(double x) const {
            return (x >= x_min && x <= x_max);
        }
        bool on_boundary(double x) const {
            return std::abs(x - x_min) < EPSILON || std::abs(x - x_max) < EPSILON;
        }
    };
    
    class Domain1d{
    private:
        Interval domain;
    public:
        bool inDomain(double x, double y) const {
            return domain.contains(x);
        }
        bool onBoundary(double x, double y) const {
            return domain.on_boundary(x);
        }
    };
}

#endif
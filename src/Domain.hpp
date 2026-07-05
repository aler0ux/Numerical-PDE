#ifndef DOMAIN_HPP
#define DOMAIN_HPP
#include <vector>
#include <variant>
#include <cmath>
namespace pde {
    constexpr double EPSILON = 1e-9;
    struct Rectangle {
        double x_min,x_max;
        double y_min,y_max;
        bool contains(double x, double y) const {
            return (x >= x_min && x <= x_max && y >= y_min && y <= y_max);
        }
        bool on_boundary(double x, double y) const {
            return contains(x,y) && (std::abs(x - x_min) < EPSILON || std::abs(x - x_max) < EPSILON || std::abs(y - y_min) < EPSILON || std::abs(y - y_max) < EPSILON);
        }
    };
    struct Circle {
        double x_center,y_center;
        double radius;
        bool contains(double x, double y) const {
            return ((x-x_center)*(x-x_center) + (y-y_center)*(y-y_center) <= radius);
        }
        bool on_boundary(double x, double y) const {
            return ((x-x_center)*(x-x_center) + (y-y_center)*(y-y_center) == radius);
        }
    };
    //TODO:add more shapes must be included in variant
    using Shape = std::variant<Rectangle,Circle>;
    enum class OpType {
        Add,
        Substract
    };
    
    struct DomainOp {
        OpType type;
        Shape shape;
    };
    class Domain{
    private:
        std::vector<DomainOp> domainOperations;
    public:
        void addShape(Shape s){
            domainOperations.push_back({OpType::Add,s});
        }
        void removeShape(Shape s){
            domainOperations.push_back({OpType::Substract,s});
        }
        bool inDomain(double x, double y) const;
        bool onBoundary(double x, double y) const;
    };
}

#endif
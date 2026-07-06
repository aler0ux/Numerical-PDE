#include "Domain2d.hpp"

namespace pde {
    template<class... Ts> struct overloaded : Ts... { using Ts::operator()...; };
    template<class... Ts> overloaded(Ts...) -> overloaded<Ts...>;
   
    bool Domain2d::inDomain(double x, double y) const {
        bool is_inside = false;
        for (const auto& op : domainOperations) {
            bool inOperation = std::visit(overloaded{
                [x,y](const Rectangle& r) { return r.contains(x,y);},
                [x,y](const Circle& r) { return r.contains(x,y);}
            }, op.shape);

            if (op.type == OpType::Add){
                is_inside = is_inside || inOperation;
            }
            else if (op.type == OpType::Substract){
                is_inside = is_inside && inOperation;
            }
        }
    }
    bool Domain2d::onBoundary(double x, double y) const {
        bool onBoundary = false;
        for (const auto& op : domainOperations) {
            bool inOperation = std::visit(overloaded{
                [x,y](const Rectangle& r) { return r.contains(x,y);},
                [x,y](const Circle& r) { return r.contains(x,y);}
            }, op.shape);
            bool onOperationBoundary = std::visit(overloaded{
                [x,y](const Rectangle& r) { return r.on_boundary(x,y);},
                [x,y](const Circle& r) { return r.on_boundary(x,y);}
            }, op.shape);
            //We now react to the matrix of possibilities

            //Point is strictly contained in an area 
            if (!onBoundary && inOperation) {
                onBoundary = false;
            }
            else if (onBoundary) { //obviously in operation also
                onBoundary = true;
            }
            //not in boundary not in operation tells us nothing about if x,y on boundary
            //in boundary not in operation makes no sense
        }
        return onBoundary;
    }


}

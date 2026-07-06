#include <iostream>
#include <functional>
#include "Expression.hpp"

namespace pde {
    //struct that defines 
    // A(du/dx) + B(du/dy) + Cu = R(x,y)
    struct Generalized2D {
    Expr A; //du / dx
    Expr B; //du / dy
    Expr C; //u
    std::function<double(double, double)> R; //function term
    };
}
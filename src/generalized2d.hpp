#include <iostream>
#include <functional>
#include "Expression.hpp"

namespace pde {
    //struct that defines 
    // A(d^2u/dx^2) + B(d^2u/dxdy) + C(d^2u/dy^2) + D(du/dx) + E(du/dy) + Fu = R(x,y)
    struct Generalized2D {
    Expr A; //d^2u / dx^2
    Expr B; //d^2u / dxdy
    Expr C; //d^2u / dy^2
    Expr D; //du   / dx
    Expr E; //du   / dy
    Expr F; //u
    std::function<double(double, double)> R; //function term
    }
}
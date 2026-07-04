#include <iostream>
#include <variant>
#include <memory>
#include <string>

namespace pde {
    //Types of sub expressions
    struct Number {double val;};
    struct Variable { std::string name;};

    //Operation declarations
    struct Add;
    struct Multiply;

    // Declare the shape of an expression
    using Expr = std::variant<Number,Variable, std::unique_ptr<Add>,std::unique_ptr<Multiply>>;

    struct Add {
        Expr left;
        Expr right;
    };
    struct Multiply{
        Expr left;
        Expr right;
    };
    

}
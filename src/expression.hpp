#ifndef EXPRESSION_HPP
#define EXPRESSION_HPP

#include <variant>
#include <memory>
#include <string>

namespace pde {
    struct Add;
    struct Multiply;

    struct Number { double val; };
    struct Variable { std::string name; };
    using Expr = std::variant<Number, Variable, std::unique_ptr<Add>, std::unique_ptr<Multiply>>;
    struct Add {
        Expr left;
        Expr right;
    };
    struct Multiply {
        Expr left;
        Expr right;
    };

    

}

#endif // EXPRESSION_HPP
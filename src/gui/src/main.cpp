#include <FL/Fl.H>

#include <memory>

#include "AppWindow.hpp"
#include "TemplateProblems.hpp"

int main(int argc, char** argv) {
    Fl::scheme("gtk+");
    Fl::visual(FL_DOUBLE | FL_RGB);

    pde::AppWindow window(1200, 760, "PDE Studio");

    // Swap these for your own IPDEProblem implementations -- the GUI does
    // not know or care what math runs inside step().
    window.addProblem(std::make_unique<pde::examples::Template1D>());
    window.addProblem(std::make_unique<pde::examples::Template2D>());

    window.show(argc, argv);
    return Fl::run();
}

#include <FL/Fl.H>

#include "CompareWindow.hpp"

int main(int argc, char** argv) {
    Fl::scheme("gtk+");
    Fl::visual(FL_DOUBLE | FL_RGB);

    cmpgui::CompareWindow window(1100, 680, "FEM vs FDM -- Poisson Comparison");
    window.show(argc, argv);
    return Fl::run();
}

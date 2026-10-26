#ifndef PDE_COMPARE_GUI_COMPAREWINDOW_HPP
#define PDE_COMPARE_GUI_COMPAREWINDOW_HPP

// The standalone interactive window: pick a domain, a source function f,
// a mesh resolution N, and a method (FDM, FEM, or both side by side),
// then press "Solve" to re-assemble and re-run Conjugate Gradient and
// refresh the heatmap(s) on the right. Nothing runs automatically when a
// control changes -- the explicit Solve button is the "refresh" step.

#include <FL/Fl_Box.H>
#include <FL/Fl_Button.H>
#include <FL/Fl_Choice.H>
#include <FL/Fl_Double_Window.H>
#include <FL/Fl_Group.H>
#include <FL/Fl_Value_Slider.H>

#include <memory>
#include <vector>

#include "../Grid2d.hpp"
#include "HeatmapView.hpp"

namespace cmpgui {

class CompareWindow : public Fl_Double_Window {
public:
    CompareWindow(int w, int h, const char* title);

private:
    static constexpr int kSidebarWidth = 300;

    Fl_Choice* domainChoice_ = nullptr;
    Fl_Choice* sourceChoice_ = nullptr;
    Fl_Choice* methodChoice_ = nullptr;
    Fl_Choice* viewChoice_ = nullptr;
    Fl_Value_Slider* nSlider_ = nullptr;
    Fl_Button* solveButton_ = nullptr;
    Fl_Box* statusBox_ = nullptr;
    HeatmapView* view_ = nullptr;

    // Kept alive between solves so HeatmapView's Panels (which only hold
    // pointers) stay valid until the next solve replaces them.
    std::unique_ptr<cmp::Grid2d> grid_;
    std::vector<double> fdmSolution_;
    std::vector<double> femSolution_;

    void buildSidebar(int x, int y, int w, int h);
    void solve();

    static void solveCB(Fl_Widget*, void* data);
    static void viewChoiceCB(Fl_Widget*, void* data);
};

}  // namespace cmpgui

#endif  // PDE_COMPARE_GUI_COMPAREWINDOW_HPP

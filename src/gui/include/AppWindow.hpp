#pragma once

#include <FL/Fl_Box.H>
#include <FL/Fl_Button.H>
#include <FL/Fl_Choice.H>
#include <FL/Fl_Double_Window.H>
#include <FL/Fl_Group.H>
#include <FL/Fl_Value_Slider.H>

#include <array>
#include <memory>
#include <string>
#include <vector>

#include "IPDEProblem.hpp"
#include "PlotWidget.hpp"

namespace pde {

// Top-level window: a sidebar of run controls + dynamically generated
// parameter sliders on the left, and the PlotWidget viewer on the right.
// Register any number of IPDEProblem implementations with addProblem();
// the GUI itself never implements or touches the numerics.
class AppWindow : public Fl_Double_Window {
public:
    AppWindow(int w, int h, const char* title);

    void addProblem(std::unique_ptr<IPDEProblem> problem);

private:
    static constexpr int kMaxParamSliders = 8;
    static constexpr int kSidebarWidth = 260;
    static constexpr int kStatusHeight = 28;

    std::vector<std::unique_ptr<IPDEProblem>> problems_;
    IPDEProblem* current_ = nullptr;
    bool playing_ = false;

    PlotWidget* plot_ = nullptr;
    Fl_Choice* problemChoice_ = nullptr;
    Fl_Button* playButton_ = nullptr;
    Fl_Value_Slider* dtSlider_ = nullptr;
    Fl_Value_Slider* speedSlider_ = nullptr;
    Fl_Group* paramGroup_ = nullptr;
    std::array<Fl_Value_Slider*, kMaxParamSliders> paramSliders_{};
    std::array<std::string, kMaxParamSliders> paramLabels_{};
    Fl_Box* statusBox_ = nullptr;

    void buildSidebar(int x, int y, int w, int h);
    void selectProblem(int index);
    void rebuildParamSliders();
    void updateStatus();
    void doStep();

    static void timerCB(void* data);
    void tick();

    static void playPauseCB(Fl_Widget*, void* data);
    static void stepCB(Fl_Widget*, void* data);
    static void resetCB(Fl_Widget*, void* data);
    static void problemChoiceCB(Fl_Widget*, void* data);
    static void paramSliderCB(Fl_Widget* w, void* data);
};

} // namespace pde

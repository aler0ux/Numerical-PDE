#pragma once

#include <FL/Fl_Widget.H>

#include "IPDEProblem.hpp"

namespace pde {

// Renders whatever IPDEProblem::kind() reports: a line plot for Scalar1D,
// or a heatmap with a colorbar for Scalar2D. Purely a viewer -- it never
// advances the simulation itself.
class PlotWidget : public Fl_Widget {
public:
    PlotWidget(int x, int y, int w, int h);

    void setProblem(IPDEProblem* problem);
    void draw() override;

private:
    IPDEProblem* problem_ = nullptr;

    void drawBackground();
    void drawLinePlot(int px, int py, int pw, int ph);
    void drawHeatmap(int px, int py, int pw, int ph);
    void drawColorbar(int cx, int cy, int cw, int ch, double lo, double hi);
    void drawFrameAndTitle(int px, int py, int pw, int ph);
    void computeRange(double& lo, double& hi) const;
};

} // namespace pde

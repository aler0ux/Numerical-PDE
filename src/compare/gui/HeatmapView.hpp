#ifndef PDE_COMPARE_GUI_HEATMAPVIEW_HPP
#define PDE_COMPARE_GUI_HEATMAPVIEW_HPP

// A small Fl_Widget that draws one or two node-value fields (FDM and/or
// FEM) side by side, each as either a flat heatmap or a shaded 3D-looking
// surface (a software axonometric projection -- no OpenGL dependency, in
// keeping with the rest of this module staying dependency-free). Pure
// viewer: it never touches the solver, it just draws whatever Panels it
// was last given, in whichever ViewMode is currently selected.

#include <FL/Fl_Widget.H>

#include <string>
#include <vector>

#include "../Grid2d.hpp"

namespace cmpgui {

struct Panel {
    std::string title;
    const cmp::Grid2d* grid = nullptr;
    const std::vector<double>* values = nullptr;  // indexed by grid DOF
};

enum class ViewMode { Heatmap, Surface3D };

class HeatmapView : public Fl_Widget {
public:
    HeatmapView(int x, int y, int w, int h);

    void setPanels(std::vector<Panel> panels);
    void setViewMode(ViewMode mode);
    ViewMode viewMode() const { return viewMode_; }

    void draw() override;

private:
    std::vector<Panel> panels_;
    ViewMode viewMode_ = ViewMode::Heatmap;

    void drawPanelHeatmap(const Panel& panel, int px, int py, int pw, int ph) const;
    void drawPanelSurface3D(const Panel& panel, int px, int py, int pw, int ph) const;
    void drawColorbar(int cx, int cy, int cw, int ch, double lo, double hi) const;
    void computeRange(const Panel& panel, double& lo, double& hi) const;
};

}  // namespace cmpgui

#endif  // PDE_COMPARE_GUI_HEATMAPVIEW_HPP

#include "HeatmapView.hpp"

#include <FL/Fl.H>
#include <FL/fl_draw.H>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <limits>

#include "../../gui/include/ColorMap.hpp"  // pde::viridis, pde::RGB -- reused as-is

namespace cmpgui {

namespace {
const Fl_Color kBackground = fl_rgb_color(24, 26, 32);
const Fl_Color kPanelBg = fl_rgb_color(32, 35, 43);
const Fl_Color kGridLine = fl_rgb_color(52, 56, 66);
const Fl_Color kMeshEdge = fl_rgb_color(18, 19, 24);
const Fl_Color kAxis = fl_rgb_color(170, 174, 184);
const Fl_Color kText = fl_rgb_color(235, 236, 240);

// Axonometric (isometric-style) projection constants for the 3D surface
// view: world (x, y, z) with x,y in [0,1] (domain coords) and z in [0,1]
// (the field value, normalised per panel) maps to a 2D point via a fixed
// 30-degree oblique projection -- the same simple technique used for
// classic "2.5D" game/engineering-drawing views, not a full 3D camera.
constexpr double kCosA = 0.8660254037844386;  // cos(30 deg)
constexpr double kSinA = 0.5;                 // sin(30 deg)
constexpr double kElevation = 0.6;            // vertical exaggeration of the (normalised) field

struct Point2 {
    double x, y;
};

Point2 project(double gx, double gy, double zNorm) {
    return Point2{(gx - gy) * kCosA, (gx + gy) * kSinA - zNorm * kElevation};
}

}  // namespace

HeatmapView::HeatmapView(int x, int y, int w, int h) : Fl_Widget(x, y, w, h) {}

void HeatmapView::setPanels(std::vector<Panel> panels) {
    panels_ = std::move(panels);
    redraw();
}

void HeatmapView::setViewMode(ViewMode mode) {
    viewMode_ = mode;
    redraw();
}

void HeatmapView::computeRange(const Panel& panel, double& lo, double& hi) const {
    lo = std::numeric_limits<double>::max();
    hi = std::numeric_limits<double>::lowest();
    if (!panel.grid || !panel.values) {
        lo = 0.0;
        hi = 1.0;
        return;
    }
    const cmp::Grid2d& grid = *panel.grid;
    for (int i = 0; i <= grid.n(); ++i) {
        for (int j = 0; j <= grid.n(); ++j) {
            if (!grid.active(i, j)) continue;
            const double v = (*panel.values)[static_cast<size_t>(grid.dof(i, j))];
            lo = std::min(lo, v);
            hi = std::max(hi, v);
        }
    }
    if (hi <= lo) hi = lo + 1.0;
}

void HeatmapView::draw() {
    fl_color(kBackground);
    fl_rectf(x(), y(), w(), h());

    if (panels_.empty()) {
        fl_color(kAxis);
        fl_font(FL_HELVETICA, 14);
        fl_draw("Choose settings on the left, then press Solve", x(), y(), w(), h(), FL_ALIGN_CENTER);
        return;
    }

    const int n = static_cast<int>(panels_.size());
    const int gap = 28;
    const int slotW = (w() - gap * (n + 1)) / n;
    const int topMargin = 34;
    const int bottomMargin = 16;
    const int rightPad = (viewMode_ == ViewMode::Heatmap) ? 54 : 10;

    for (int k = 0; k < n; ++k) {
        const int px = x() + gap + k * (slotW + gap);
        const int py = y() + topMargin;
        const int pw = std::max(10, slotW - rightPad);
        const int ph = std::max(10, h() - topMargin - bottomMargin);
        const Panel& panel = panels_[static_cast<size_t>(k)];

        fl_font(FL_HELVETICA_BOLD, 13);
        fl_color(kText);
        fl_draw(panel.title.c_str(), px, py - 20, pw, 18, FL_ALIGN_LEFT);

        if (viewMode_ == ViewMode::Heatmap) {
            drawPanelHeatmap(panel, px, py, pw, ph);
        } else {
            drawPanelSurface3D(panel, px, py, pw, ph);
        }
    }
}

void HeatmapView::drawPanelHeatmap(const Panel& panel, int px, int py, int pw, int ph) const {
    fl_color(kPanelBg);
    fl_rectf(px, py, pw, ph);

    if (!panel.grid || !panel.values) {
        fl_color(kGridLine);
        fl_rect(px, py, pw, ph);
        return;
    }

    const cmp::Grid2d& grid = *panel.grid;
    const int side = grid.side();

    double lo, hi;
    computeRange(panel, lo, hi);

    const double cellW = static_cast<double>(pw) / side;
    const double cellH = static_cast<double>(ph) / side;

    for (int j = 0; j < side; ++j) {
        // Flip vertically: the domain's y=0 sits at the bottom of the panel,
        // matching the plot_compare.py convention, not raw row order.
        const int row = side - 1 - j;
        const int cy = py + static_cast<int>(row * cellH);
        const int ch = std::max(1, static_cast<int>((row + 1) * cellH) - static_cast<int>(row * cellH));
        for (int i = 0; i < side; ++i) {
            if (!grid.active(i, j)) continue;
            const double v = (*panel.values)[static_cast<size_t>(grid.dof(i, j))];
            const double t = (v - lo) / (hi - lo);
            const pde::RGB c = pde::viridis(t);

            const int cx = px + static_cast<int>(i * cellW);
            const int cw = std::max(1, static_cast<int>((i + 1) * cellW) - static_cast<int>(i * cellW));
            fl_color(fl_rgb_color(c.r, c.g, c.b));
            fl_rectf(cx, cy, cw, ch);
        }
    }

    fl_color(kGridLine);
    fl_rect(px, py, pw, ph);

    drawColorbar(px + pw + 10, py, 16, ph, lo, hi);
}

void HeatmapView::drawPanelSurface3D(const Panel& panel, int px, int py, int pw, int ph) const {
    fl_color(kPanelBg);
    fl_rectf(px, py, pw, ph);
    fl_color(kGridLine);
    fl_rect(px, py, pw, ph);

    if (!panel.grid || !panel.values) return;
    const cmp::Grid2d& grid = *panel.grid;
    const int n = grid.n();

    double lo, hi;
    computeRange(panel, lo, hi);

    // Fit the projected unit cube (x,y in [0,1], z in [0,1]) into the panel.
    double minX = std::numeric_limits<double>::max(), maxX = std::numeric_limits<double>::lowest();
    double minY = std::numeric_limits<double>::max(), maxY = std::numeric_limits<double>::lowest();
    for (double gx : {0.0, 1.0})
        for (double gy : {0.0, 1.0})
            for (double zn : {0.0, 1.0}) {
                const Point2 p = project(gx, gy, zn);
                minX = std::min(minX, p.x);
                maxX = std::max(maxX, p.x);
                minY = std::min(minY, p.y);
                maxY = std::max(maxY, p.y);
            }
    const double spanX = std::max(1e-9, maxX - minX);
    const double spanY = std::max(1e-9, maxY - minY);
    const double scale = 0.92 * std::min(pw / spanX, ph / spanY);
    const double marginX = (pw - spanX * scale) / 2.0;
    const double marginY = (ph - spanY * scale) / 2.0;

    auto toScreen = [&](double gx, double gy, double zn) -> std::array<double, 2> {
        const Point2 p = project(gx, gy, zn);
        return {px + marginX + (p.x - minX) * scale, py + marginY + (p.y - minY) * scale};
    };

    auto heightAt = [&](int i, int j) -> double {
        if (!grid.active(i, j)) return 0.0;
        const double v = (*panel.values)[static_cast<size_t>(grid.dof(i, j))];
        return (v - lo) / (hi - lo);
    };

    // Painter's algorithm: draw cells back-to-front. With this projection,
    // depth increases with i+j, so drawing in that order is a correct (for
    // a reasonably smooth, single-blob field) approximation without
    // needing a real z-buffer -- consistent with the "simple technique"
    // used throughout this module.
    for (int depth = 0; depth <= 2 * (n - 1); ++depth) {
        for (int i = 0; i < n; ++i) {
            const int j = depth - i;
            if (j < 0 || j >= n) continue;
            const bool bl = grid.active(i, j);
            const bool br = grid.active(i + 1, j);
            const bool tr = grid.active(i + 1, j + 1);
            const bool tl = grid.active(i, j + 1);
            if (!(bl && br && tr && tl)) continue;  // only fully-active cells, same as the FEM mesh's 4-corner case

            const double x0 = grid.x(i), x1 = grid.x(i + 1);
            const double y0 = grid.y(j), y1 = grid.y(j + 1);
            const double zbl = heightAt(i, j), zbr = heightAt(i + 1, j);
            const double ztr = heightAt(i + 1, j + 1), ztl = heightAt(i, j + 1);

            const auto sBl = toScreen(x0, y0, zbl);
            const auto sBr = toScreen(x1, y0, zbr);
            const auto sTr = toScreen(x1, y1, ztr);
            const auto sTl = toScreen(x0, y1, ztl);

            const double avgZ = (zbl + zbr + ztr + ztl) / 4.0;
            const pde::RGB c = pde::viridis(avgZ);
            fl_color(fl_rgb_color(c.r, c.g, c.b));
            fl_begin_polygon();
            fl_vertex(sBl[0], sBl[1]);
            fl_vertex(sBr[0], sBr[1]);
            fl_vertex(sTr[0], sTr[1]);
            fl_vertex(sTl[0], sTl[1]);
            fl_end_polygon();

            fl_color(kMeshEdge);
            fl_begin_loop();
            fl_vertex(sBl[0], sBl[1]);
            fl_vertex(sBr[0], sBr[1]);
            fl_vertex(sTr[0], sTr[1]);
            fl_vertex(sTl[0], sTl[1]);
            fl_end_loop();
        }
    }
}

void HeatmapView::drawColorbar(int cx, int cy, int cw, int ch, double lo, double hi) const {
    const int safeCh = std::max(ch, 2);
    for (int row = 0; row < safeCh; ++row) {
        const double t = 1.0 - static_cast<double>(row) / (safeCh - 1);
        const pde::RGB c = pde::viridis(t);
        fl_color(fl_rgb_color(c.r, c.g, c.b));
        fl_line(cx, cy + row, cx + cw, cy + row);
    }
    fl_color(kGridLine);
    fl_rect(cx, cy, cw, ch);

    fl_font(FL_HELVETICA, 10);
    fl_color(kAxis);
    char buf[32];
    std::snprintf(buf, sizeof(buf), "%.3g", hi);
    fl_draw(buf, cx - 4, cy - 2, cw + 44, 12, FL_ALIGN_LEFT);
    std::snprintf(buf, sizeof(buf), "%.3g", lo);
    fl_draw(buf, cx - 4, cy + ch - 10, cw + 44, 12, FL_ALIGN_LEFT);
}

}  // namespace cmpgui

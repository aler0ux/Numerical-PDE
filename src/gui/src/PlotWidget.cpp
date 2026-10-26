#include "PlotWidget.hpp"

#include <FL/Fl.H>
#include <FL/fl_draw.H>

#include <algorithm>
#include <cstdio>
#include <limits>

#include "ColorMap.hpp"

namespace pde {

namespace {
const Fl_Color kBackground = fl_rgb_color(24, 26, 32);
const Fl_Color kPanel = fl_rgb_color(32, 35, 43);
const Fl_Color kGrid = fl_rgb_color(52, 56, 66);
const Fl_Color kAxis = fl_rgb_color(140, 146, 158);
const Fl_Color kText = fl_rgb_color(220, 222, 228);
const Fl_Color kCurve = fl_rgb_color(96, 200, 255);
} // namespace

PlotWidget::PlotWidget(int x, int y, int w, int h) : Fl_Widget(x, y, w, h) {}

void PlotWidget::setProblem(IPDEProblem* problem) {
    problem_ = problem;
    redraw();
}

void PlotWidget::computeRange(double& lo, double& hi) const {
    lo = std::numeric_limits<double>::max();
    hi = std::numeric_limits<double>::lowest();
    if (!problem_) {
        lo = 0.0;
        hi = 1.0;
        return;
    }

    if (double flo, fhi; problem_->fixedRange(flo, fhi)) {
        lo = flo;
        hi = fhi;
        return;
    }

    const int n = problem_->nx();
    const int m = problem_->kind() == FieldKind::Scalar2D ? problem_->ny() : 1;
    for (int j = 0; j < m; ++j) {
        for (int i = 0; i < n; ++i) {
            const double v = problem_->valueAt(i, j);
            lo = std::min(lo, v);
            hi = std::max(hi, v);
        }
    }
    if (hi <= lo) {
        hi = lo + 1.0;
    }
}

void PlotWidget::drawBackground() {
    fl_color(kBackground);
    fl_rectf(x(), y(), w(), h());
}

void PlotWidget::draw() {
    drawBackground();

    if (!problem_) {
        fl_color(kAxis);
        fl_font(FL_HELVETICA, 14);
        fl_draw("No problem selected", x(), y(), w(), h(), FL_ALIGN_CENTER);
        return;
    }

    const bool is2D = problem_->kind() == FieldKind::Scalar2D;

    const int leftMargin = 56;
    const int rightMargin = is2D ? 84 : 24;
    const int topMargin = 40;
    const int bottomMargin = 40;

    const int px = x() + leftMargin;
    const int py = y() + topMargin;
    const int pw = std::max(10, w() - leftMargin - rightMargin);
    const int ph = std::max(10, h() - topMargin - bottomMargin);

    fl_color(kPanel);
    fl_rectf(px, py, pw, ph);

    if (is2D) {
        drawHeatmap(px, py, pw, ph);
    } else {
        drawLinePlot(px, py, pw, ph);
    }

    drawFrameAndTitle(px, py, pw, ph);
}

void PlotWidget::drawFrameAndTitle(int px, int py, int pw, int ph) {
    fl_color(kGrid);
    fl_rect(px, py, pw, ph);

    fl_color(kText);
    fl_font(FL_HELVETICA_BOLD, 15);
    fl_draw(problem_->name().c_str(), x() + 12, y() + 24);

    char buf[96];
    std::snprintf(buf, sizeof(buf), "t = %.4f", problem_->time());
    fl_font(FL_HELVETICA, 12);
    fl_color(kAxis);
    fl_draw(buf, px, y() + 24, pw, 20, FL_ALIGN_RIGHT);
}

void PlotWidget::drawLinePlot(int px, int py, int pw, int ph) {
    const int n = problem_->nx();
    if (n < 2) {
        return;
    }

    double lo, hi;
    computeRange(lo, hi);

    // Horizontal gridlines + value labels.
    fl_font(FL_HELVETICA, 10);
    for (int k = 0; k <= 4; ++k) {
        const double frac = static_cast<double>(k) / 4.0;
        const int gy = py + static_cast<int>(frac * ph);
        fl_color(kGrid);
        fl_line(px, gy, px + pw, gy);

        char buf[32];
        std::snprintf(buf, sizeof(buf), "%.3g", hi - frac * (hi - lo));
        fl_color(kAxis);
        fl_draw(buf, px - 52, gy - 6, 48, 12, FL_ALIGN_RIGHT);
    }

    // The field as a connected polyline.
    fl_color(kCurve);
    int prevX = 0, prevY = 0;
    for (int i = 0; i < n; ++i) {
        const double v = problem_->valueAt(i);
        const int sx = px + static_cast<int>((static_cast<double>(i) / (n - 1)) * pw);
        const int sy = py + ph - static_cast<int>(((v - lo) / (hi - lo)) * ph);
        if (i > 0) {
            fl_line(prevX, prevY, sx, sy);
        }
        prevX = sx;
        prevY = sy;
    }
}

void PlotWidget::drawHeatmap(int px, int py, int pw, int ph) {
    const int n = problem_->nx();
    const int m = problem_->ny();
    if (n < 1 || m < 1) {
        return;
    }

    double lo, hi;
    computeRange(lo, hi);

    const double cellW = static_cast<double>(pw) / n;
    const double cellH = static_cast<double>(ph) / m;

    for (int j = 0; j < m; ++j) {
        const int cy = py + static_cast<int>(j * cellH);
        const int ch = std::max(1, static_cast<int>((j + 1) * cellH) - static_cast<int>(j * cellH));
        for (int i = 0; i < n; ++i) {
            const double v = problem_->valueAt(i, j);
            const double t = (v - lo) / (hi - lo);
            const RGB c = viridis(t);

            const int cx = px + static_cast<int>(i * cellW);
            const int cw = std::max(1, static_cast<int>((i + 1) * cellW) - static_cast<int>(i * cellW));
            fl_color(fl_rgb_color(c.r, c.g, c.b));
            fl_rectf(cx, cy, cw, ch);
        }
    }

    drawColorbar(px + pw + 16, py, 20, ph, lo, hi);
}

void PlotWidget::drawColorbar(int cx, int cy, int cw, int ch, double lo, double hi) {
    for (int row = 0; row < ch; ++row) {
        const double t = 1.0 - static_cast<double>(row) / (ch - 1);
        const RGB c = viridis(t);
        fl_color(fl_rgb_color(c.r, c.g, c.b));
        fl_line(cx, cy + row, cx + cw, cy + row);
    }
    fl_color(kGrid);
    fl_rect(cx, cy, cw, ch);

    fl_font(FL_HELVETICA, 10);
    fl_color(kAxis);
    char buf[32];
    std::snprintf(buf, sizeof(buf), "%.3g", hi);
    fl_draw(buf, cx - 4, cy - 2, cw + 40, 12, FL_ALIGN_LEFT);
    std::snprintf(buf, sizeof(buf), "%.3g", lo);
    fl_draw(buf, cx - 4, cy + ch - 10, cw + 40, 12, FL_ALIGN_LEFT);
}

} // namespace pde

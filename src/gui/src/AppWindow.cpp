#include "AppWindow.hpp"

#include <FL/Fl.H>
#include <FL/fl_draw.H>

#include <cstdio>
#include <limits>

namespace pde {

namespace {
const Fl_Color kWindowBg = fl_rgb_color(20, 22, 27);
const Fl_Color kSidebarBg = fl_rgb_color(28, 30, 37);
const Fl_Color kText = fl_rgb_color(220, 222, 228);
const Fl_Color kMuted = fl_rgb_color(150, 155, 165);
const Fl_Color kAccent = fl_rgb_color(96, 200, 255);
constexpr double kFrameInterval = 1.0 / 30.0;
} // namespace

AppWindow::AppWindow(int w, int h, const char* title) : Fl_Double_Window(w, h, title) {
    color(kWindowBg);

    buildSidebar(0, 0, kSidebarWidth, h - kStatusHeight);

    plot_ = new PlotWidget(kSidebarWidth, 0, w - kSidebarWidth, h - kStatusHeight);

    statusBox_ = new Fl_Box(0, h - kStatusHeight, w, kStatusHeight, "");
    statusBox_->box(FL_FLAT_BOX);
    statusBox_->color(kSidebarBg);
    statusBox_->labelcolor(kMuted);
    statusBox_->labelsize(12);
    statusBox_->align(FL_ALIGN_LEFT | FL_ALIGN_INSIDE);

    end();
    resizable(plot_);
}

void AppWindow::buildSidebar(int x, int y, int w, int h) {
    Fl_Group* sidebar = new Fl_Group(x, y, w, h);
    sidebar->box(FL_FLAT_BOX);
    sidebar->color(kSidebarBg);

    const int pad = 14;
    int cy = y + pad;
    const int cw = w - 2 * pad;

    auto* title = new Fl_Box(x + pad, cy, cw, 26, "PDE Studio");
    title->labelfont(FL_HELVETICA_BOLD);
    title->labelsize(18);
    title->labelcolor(kText);
    title->align(FL_ALIGN_LEFT | FL_ALIGN_INSIDE);
    cy += 34;

    problemChoice_ = new Fl_Choice(x + pad, cy, cw, 26);
    problemChoice_->callback(problemChoiceCB, this);
    cy += 38;

    playButton_ = new Fl_Button(x + pad, cy, cw / 2 - 4, 30, "Play");
    playButton_->callback(playPauseCB, this);
    auto* stepButton = new Fl_Button(x + pad + cw / 2 + 4, cy, cw / 2 - 4, 30, "Step");
    stepButton->callback(stepCB, this);
    cy += 38;

    auto* resetButton = new Fl_Button(x + pad, cy, cw, 26, "Reset");
    resetButton->callback(resetCB, this);
    cy += 40;

    auto* dtLabel = new Fl_Box(x + pad, cy, cw, 16, "Time step (dt)");
    dtLabel->labelsize(11);
    dtLabel->labelcolor(kMuted);
    dtLabel->align(FL_ALIGN_LEFT | FL_ALIGN_INSIDE);
    cy += 18;
    dtSlider_ = new Fl_Value_Slider(x + pad, cy, cw, 22);
    dtSlider_->type(FL_HOR_NICE_SLIDER);
    dtSlider_->bounds(0.0005, 0.05);
    dtSlider_->value(0.005);
    dtSlider_->step(0.0005);
    cy += 32;

    auto* speedLabel = new Fl_Box(x + pad, cy, cw, 16, "Steps per frame");
    speedLabel->labelsize(11);
    speedLabel->labelcolor(kMuted);
    speedLabel->align(FL_ALIGN_LEFT | FL_ALIGN_INSIDE);
    cy += 18;
    speedSlider_ = new Fl_Value_Slider(x + pad, cy, cw, 22);
    speedSlider_->type(FL_HOR_NICE_SLIDER);
    speedSlider_->bounds(1, 20);
    speedSlider_->value(1);
    speedSlider_->step(1);
    cy += 34;

    auto* sep = new Fl_Box(x + pad, cy, cw, 1, "");
    sep->box(FL_FLAT_BOX);
    sep->color(fl_rgb_color(60, 64, 74));
    cy += 10;

    auto* paramTitle = new Fl_Box(x + pad, cy, cw, 18, "Parameters");
    paramTitle->labelfont(FL_HELVETICA_BOLD);
    paramTitle->labelsize(12);
    paramTitle->labelcolor(kText);
    paramTitle->align(FL_ALIGN_LEFT | FL_ALIGN_INSIDE);
    cy += 24;

    paramGroup_ = new Fl_Group(x + pad, cy, cw, h - (cy - y) - pad);
    int py = cy;
    for (int i = 0; i < kMaxParamSliders; ++i) {
        auto* slider = new Fl_Value_Slider(x + pad, py, cw, 22);
        slider->type(FL_HOR_NICE_SLIDER);
        slider->callback(paramSliderCB, this);
        slider->hide();
        paramSliders_[i] = slider;
        py += 30;
    }
    paramGroup_->end();

    sidebar->end();

    for (int i = 0; i < sidebar->children(); ++i) {
        Fl_Widget* child = sidebar->child(i);
        if (dynamic_cast<Fl_Button*>(child) || dynamic_cast<Fl_Choice*>(child) ||
            dynamic_cast<Fl_Value_Slider*>(child)) {
            child->labelcolor(kText);
            child->color(fl_rgb_color(46, 49, 58));
            child->selection_color(kAccent);
        }
    }
    // Sidebar's own controls also need the same styling; children() above
    // only covers direct children, so restyle the nested param sliders too.
    for (int i = 0; i < kMaxParamSliders; ++i) {
        paramSliders_[i]->color(fl_rgb_color(46, 49, 58));
        paramSliders_[i]->selection_color(kAccent);
        paramSliders_[i]->labelcolor(kText);
    }
}

void AppWindow::addProblem(std::unique_ptr<IPDEProblem> problem) {
    problemChoice_->add(problem->name().c_str());
    problems_.push_back(std::move(problem));
    if (problems_.size() == 1) {
        problemChoice_->value(0);
        selectProblem(0);
    }
}

void AppWindow::selectProblem(int index) {
    if (index < 0 || static_cast<size_t>(index) >= problems_.size()) {
        return;
    }
    current_ = problems_[static_cast<size_t>(index)].get();
    current_->reset();
    plot_->setProblem(current_);
    rebuildParamSliders();
    updateStatus();
    plot_->redraw();
}

void AppWindow::rebuildParamSliders() {
    auto& params = current_->parameters();
    for (int i = 0; i < kMaxParamSliders; ++i) {
        if (static_cast<size_t>(i) < params.size()) {
            const Parameter& p = params[static_cast<size_t>(i)];
            paramLabels_[i] = p.name;
            Fl_Value_Slider* s = paramSliders_[i];
            s->bounds(p.min, p.max);
            s->step(p.step > 0.0 ? p.step : (p.max - p.min) / 200.0);
            s->value(p.value);
            s->copy_label(paramLabels_[i].c_str());
            s->align(FL_ALIGN_TOP_LEFT);
            s->user_data(reinterpret_cast<void*>(static_cast<intptr_t>(i)));
            s->show();
        } else {
            paramSliders_[i]->hide();
        }
    }
    paramGroup_->redraw();
}

void AppWindow::updateStatus() {
    if (!current_) {
        statusBox_->label("");
        return;
    }
    double lo = std::numeric_limits<double>::max();
    double hi = std::numeric_limits<double>::lowest();
    const int n = current_->nx();
    const int m = current_->kind() == FieldKind::Scalar2D ? current_->ny() : 1;
    for (int j = 0; j < m; ++j) {
        for (int i = 0; i < n; ++i) {
            const double v = current_->valueAt(i, j);
            if (v < lo) lo = v;
            if (v > hi) hi = v;
        }
    }

    const std::string gridStr = (m > 1) ? (std::to_string(n) + " x " + std::to_string(m)) : std::to_string(n);
    static char buf[160];
    std::snprintf(buf, sizeof(buf), "  t = %.4f    grid = %s    min = %.4g    max = %.4g", current_->time(),
                  gridStr.c_str(), lo, hi);
    statusBox_->copy_label(buf);
}

void AppWindow::doStep() {
    if (!current_) {
        return;
    }
    const int steps = static_cast<int>(speedSlider_->value());
    const double dt = dtSlider_->value();
    for (int i = 0; i < steps; ++i) {
        current_->step(dt);
    }
    updateStatus();
    plot_->redraw();
}

void AppWindow::tick() {
    if (playing_) {
        doStep();
        Fl::repeat_timeout(kFrameInterval, timerCB, this);
    }
}

void AppWindow::timerCB(void* data) {
    static_cast<AppWindow*>(data)->tick();
}

void AppWindow::playPauseCB(Fl_Widget* w, void* data) {
    auto* self = static_cast<AppWindow*>(data);
    self->playing_ = !self->playing_;
    static_cast<Fl_Button*>(w)->copy_label(self->playing_ ? "Pause" : "Play");
    if (self->playing_) {
        Fl::add_timeout(kFrameInterval, timerCB, self);
    }
}

void AppWindow::stepCB(Fl_Widget*, void* data) {
    static_cast<AppWindow*>(data)->doStep();
}

void AppWindow::resetCB(Fl_Widget*, void* data) {
    auto* self = static_cast<AppWindow*>(data);
    if (self->current_) {
        self->current_->reset();
        self->updateStatus();
        self->plot_->redraw();
    }
}

void AppWindow::problemChoiceCB(Fl_Widget*, void* data) {
    auto* self = static_cast<AppWindow*>(data);
    self->selectProblem(self->problemChoice_->value());
}

void AppWindow::paramSliderCB(Fl_Widget* w, void* data) {
    auto* self = static_cast<AppWindow*>(data);
    if (!self->current_) {
        return;
    }
    const auto index = static_cast<int>(reinterpret_cast<intptr_t>(w->user_data()));
    auto& params = self->current_->parameters();
    if (index >= 0 && static_cast<size_t>(index) < params.size()) {
        params[static_cast<size_t>(index)].value = static_cast<Fl_Value_Slider*>(w)->value();
    }
}

} // namespace pde

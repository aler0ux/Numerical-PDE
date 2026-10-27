#include "CompareWindow.hpp"

#include <FL/Fl.H>
#include <FL/fl_draw.H>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <sstream>

#include "../Assemble.hpp"
#include "../ConjugateGradient.hpp"
#include "../Domains.hpp"
#include "../SourceFunctions.hpp"

namespace cmpgui {

namespace {
const Fl_Color kWindowBg = fl_rgb_color(20, 22, 27);
const Fl_Color kSidebarBg = fl_rgb_color(28, 30, 37);
// Pure/near-white throughout: every label, dropdown entry, and slider
// readout sits on a dark background, so anything less than this reads as
// illegible grey-on-charcoal.
const Fl_Color kText = FL_WHITE;
const Fl_Color kMuted = fl_rgb_color(205, 208, 215);
const Fl_Color kAccent = fl_rgb_color(96, 200, 255);
const Fl_Color kControlBg = fl_rgb_color(46, 49, 58);

// Fixed display order for the domain dropdown; index into this array is
// exactly Fl_Choice::value() for domainChoice_.
constexpr cmp::DomainKind kDomainOrder[] = {cmp::DomainKind::Square, cmp::DomainKind::LShape,
                                             cmp::DomainKind::Plus, cmp::DomainKind::Circle,
                                             cmp::DomainKind::Eight, cmp::DomainKind::Donut};

enum class Method { Fdm, Fem, Both };

}  // namespace

CompareWindow::CompareWindow(int w, int h, const char* title) : Fl_Double_Window(w, h, title) {
    color(kWindowBg);

    buildSidebar(0, 0, kSidebarWidth, h);
    view_ = new HeatmapView(kSidebarWidth, 0, w - kSidebarWidth, h);

    end();
    resizable(view_);

    solve();  // populate the view on first open instead of starting blank
}

void CompareWindow::buildSidebar(int x, int y, int w, int h) {
    auto* sidebar = new Fl_Group(x, y, w, h);
    sidebar->box(FL_FLAT_BOX);
    sidebar->color(kSidebarBg);

    const int pad = 16;
    int cy = y + pad;
    const int cw = w - 2 * pad;

    auto* title = new Fl_Box(x + pad, cy, cw, 28, "FEM vs FDM: Poisson");
    title->labelfont(FL_HELVETICA_BOLD);
    title->labelsize(16);
    title->labelcolor(kText);
    title->align(FL_ALIGN_LEFT | FL_ALIGN_INSIDE);
    cy += 38;

    auto label = [&](const char* text) {
        auto* box = new Fl_Box(x + pad, cy, cw, 16, text);
        box->labelsize(11);
        box->labelcolor(kMuted);
        box->align(FL_ALIGN_LEFT | FL_ALIGN_INSIDE);
        cy += 18;
    };

    label("Domain");
    domainChoice_ = new Fl_Choice(x + pad, cy, cw, 26);
    for (cmp::DomainKind d : kDomainOrder) domainChoice_->add(cmp::domainName(d).c_str());
    domainChoice_->value(3);  // circle: the most visually interesting default
    cy += 36;

    label("Source function f(x,y)");
    sourceChoice_ = new Fl_Choice(x + pad, cy, cw, 26);
    int sineIndex = 0;
    const auto& library = cmp::sourceFunctionLibrary();
    for (size_t k = 0; k < library.size(); ++k) {
        sourceChoice_->add(library[k].name.c_str());
        if (library[k].name == "sine") sineIndex = static_cast<int>(k);
    }
    sourceChoice_->value(sineIndex);  // "sine" -- found by name, not a hardcoded
                                       // index, so it stays correct as the library grows.
    cy += 36;

    label("Method");
    methodChoice_ = new Fl_Choice(x + pad, cy, cw, 26);
    methodChoice_->add("FDM only");
    methodChoice_->add("FEM only");
    methodChoice_->add("FDM vs FEM (both)");
    methodChoice_->value(2);
    cy += 36;

    label("Mesh resolution N -- FDM (cells per axis)");
    nSliderFdm_ = new Fl_Value_Slider(x + pad, cy, cw, 24);
    nSliderFdm_->type(FL_HOR_NICE_SLIDER);
    nSliderFdm_->bounds(8, 160);
    nSliderFdm_->value(64);
    nSliderFdm_->step(1);
    cy += 36;

    label("Mesh resolution N -- FEM (cells per axis)");
    nSliderFem_ = new Fl_Value_Slider(x + pad, cy, cw, 24);
    nSliderFem_->type(FL_HOR_NICE_SLIDER);
    nSliderFem_->bounds(8, 160);
    nSliderFem_->value(64);
    nSliderFem_->step(1);
    cy += 36;

    label("View");
    viewChoice_ = new Fl_Choice(x + pad, cy, cw, 26);
    viewChoice_->add("Heatmap");
    viewChoice_->add("3D surface");
    viewChoice_->value(0);
    viewChoice_->callback(viewChoiceCB, this);
    cy += 36;

    solveButton_ = new Fl_Button(x + pad, cy, cw, 32, "Solve");
    solveButton_->callback(solveCB, this);
    cy += 44;

    auto* sep = new Fl_Box(x + pad, cy, cw, 1, "");
    sep->box(FL_FLAT_BOX);
    sep->color(fl_rgb_color(60, 64, 74));
    cy += 12;

    statusBox_ = new Fl_Box(x + pad, cy, cw, h - (cy - y) - pad, "");
    statusBox_->labelsize(12);
    statusBox_->labelcolor(kMuted);
    statusBox_->align(FL_ALIGN_LEFT | FL_ALIGN_TOP | FL_ALIGN_INSIDE | FL_ALIGN_WRAP);

    sidebar->end();

    for (int i = 0; i < sidebar->children(); ++i) {
        Fl_Widget* child = sidebar->child(i);
        if (dynamic_cast<Fl_Button*>(child)) {
            child->labelcolor(kText);
            child->color(kControlBg);
            child->selection_color(kAccent);
        } else if (auto* choice = dynamic_cast<Fl_Choice*>(child)) {
            choice->labelcolor(kText);
            choice->color(kControlBg);
            choice->selection_color(kAccent);
            choice->textcolor(kText);  // the selected item's own text, drawn via textcolor() not labelcolor()
        } else if (auto* slider = dynamic_cast<Fl_Value_Slider*>(child)) {
            slider->labelcolor(kText);
            slider->color(kControlBg);
            slider->selection_color(kAccent);
            slider->textcolor(kText);  // the numeric readout box
        }
    }
}

void CompareWindow::solve() {
    const cmp::DomainKind domain = kDomainOrder[static_cast<size_t>(std::max(0, domainChoice_->value()))];
    const cmp::SourceFunction& source =
        cmp::sourceFunctionLibrary()[static_cast<size_t>(std::max(0, sourceChoice_->value()))];
    const auto method = static_cast<Method>(std::max(0, methodChoice_->value()));
    const int nFdm = std::max(2, static_cast<int>(nSliderFdm_->value()));
    const int nFem = std::max(2, static_cast<int>(nSliderFem_->value()));

    std::vector<Panel> panels;
    std::ostringstream status;
    status << "Domain: " << cmp::domainName(domain) << "\n";
    status << "f = " << source.name << "\n\n";

    auto runMethod = [&](const char* label, const cmp::Grid2d& grid, const cmp::System& sys,
                          std::vector<double>& storeInto) -> cmp::CGResult {
        const auto t0 = std::chrono::steady_clock::now();
        cmp::CGResult result = cmp::conjugateGradient(sys.A, sys.b);
        const auto t1 = std::chrono::steady_clock::now();
        const double ms = std::chrono::duration<double, std::milli>(t1 - t0).count();
        storeInto = result.x;
        status << label << " (N=" << grid.n() << "):\n";
        status << "  active DOFs = " << grid.numDofs() << "\n";
        status << "  nonzeros   = " << sys.A.nnz() << "\n";
        status << "  CG iters   = " << result.iterations << (result.converged ? " (ok)" : " (NOT converged)") << "\n";
        status << "  resid      = " << result.residualHistory.back() << "\n";
        status << "  solve time = " << ms << " ms\n\n";
        return result;
    };

    if (method == Method::Fdm || method == Method::Both) {
        gridFdm_ = std::make_unique<cmp::Grid2d>(domain, nFdm);
        const cmp::System sys = cmp::assembleFDM(*gridFdm_, source.expr);
        runMethod("Finite differences", *gridFdm_, sys, fdmSolution_);
        panels.push_back(Panel{"FDM", gridFdm_.get(), &fdmSolution_});
    }
    if (method == Method::Fem || method == Method::Both) {
        gridFem_ = std::make_unique<cmp::Grid2d>(domain, nFem);
        const cmp::System sys = cmp::assembleFEM(*gridFem_, source.expr);
        runMethod("Finite elements", *gridFem_, sys, femSolution_);
        panels.push_back(Panel{"FEM", gridFem_.get(), &femSolution_});
    }

    if (method == Method::Both) {
        // Only a meaningful comparison when both methods share the same N --
        // that's what makes their DOF numbering (and hence node positions)
        // line up index-for-index (Grid2d::build() is a pure function of
        // (domain, n), so matching n gives matching grids). With
        // independent sliders the two N's can now differ, in which case
        // fdmSolution_[d] and femSolution_[d] are not even the same
        // physical point, so comparing them directly would be meaningless
        // rather than just imprecise -- skip it and say why instead.
        if (nFdm == nFem) {
            double maxDiff = 0.0;
            for (int d = 0; d < gridFdm_->numDofs(); ++d)
                maxDiff = std::max(maxDiff, std::fabs(fdmSolution_[static_cast<size_t>(d)] - femSolution_[static_cast<size_t>(d)]));
            status << "Max |u_fdm - u_fem| = " << maxDiff << "\n";
        } else {
            status << "(FDM/FEM comparison needs matching N -- currently " << nFdm << " vs " << nFem << ")\n";
        }
    }

    statusBox_->copy_label(status.str().c_str());
    statusBox_->redraw();
    view_->setPanels(std::move(panels));
}

void CompareWindow::solveCB(Fl_Widget*, void* data) {
    static_cast<CompareWindow*>(data)->solve();
}

void CompareWindow::viewChoiceCB(Fl_Widget*, void* data) {
    auto* self = static_cast<CompareWindow*>(data);
    // Purely a render-mode switch -- no need to re-solve, just redraw the
    // already-computed field(s) in the other mode.
    self->view_->setViewMode(self->viewChoice_->value() == 0 ? ViewMode::Heatmap : ViewMode::Surface3D);
}

}  // namespace cmpgui

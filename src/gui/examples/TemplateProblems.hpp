#pragma once

// These are intentionally empty scaffolds, not solvers. Each one only sets
// up a grid and an initial condition so the GUI has something to display
// out of the box; step() does no numerics. Copy one of these and fill in
// step() with your own scheme (finite differences, finite elements,
// spectral, whatever) to turn it into a real PDE demo.

#include <cmath>
#include <vector>

#include "IPDEProblem.hpp"

namespace pde::examples {

class Template1D : public IPDEProblem {
public:
    explicit Template1D(int n = 200) : n_(n), field_(static_cast<size_t>(n), 0.0) { reset(); }

    std::string name() const override { return "Template: 1D Field"; }
    FieldKind kind() const override { return FieldKind::Scalar1D; }
    int nx() const override { return n_; }
    double time() const override { return t_; }

    void step(double dt) override {
        // TODO: replace with your time-marching update, e.g.
        //   for (int i = 1; i < n_ - 1; ++i)
        //       next[i] = field_[i] + dt * kappa * (field_[i-1] - 2*field_[i] + field_[i+1]) / (dx*dx);
        t_ += dt;
    }

    void reset() override {
        t_ = 0.0;
        for (int i = 0; i < n_; ++i) {
            const double x = static_cast<double>(i) / (n_ - 1);
            field_[static_cast<size_t>(i)] = std::exp(-40.0 * (x - 0.5) * (x - 0.5)); // initial condition only
        }
    }

    double valueAt(int i, int /*j*/ = 0) const override { return field_[static_cast<size_t>(i)]; }
    std::vector<Parameter>& parameters() override { return params_; }

private:
    int n_;
    double t_ = 0.0;
    std::vector<double> field_;
    std::vector<Parameter> params_{{"kappa", 0.1, 0.0, 1.0, 0.01}};
};

class Template2D : public IPDEProblem {
public:
    explicit Template2D(int n = 80) : n_(n), field_(static_cast<size_t>(n) * static_cast<size_t>(n), 0.0) {
        reset();
    }

    std::string name() const override { return "Template: 2D Field"; }
    FieldKind kind() const override { return FieldKind::Scalar2D; }
    int nx() const override { return n_; }
    int ny() const override { return n_; }
    double time() const override { return t_; }

    void step(double dt) override {
        // TODO: replace with your time-marching update, e.g. a 5-point
        // Laplacian stencil for the heat equation or a leapfrog scheme
        // for the wave equation.
        t_ += dt;
    }

    void reset() override {
        t_ = 0.0;
        for (int j = 0; j < n_; ++j) {
            for (int i = 0; i < n_; ++i) {
                const double x = static_cast<double>(i) / (n_ - 1) - 0.5;
                const double y = static_cast<double>(j) / (n_ - 1) - 0.5;
                field_[static_cast<size_t>(j) * static_cast<size_t>(n_) + static_cast<size_t>(i)] =
                    std::exp(-40.0 * (x * x + y * y)); // initial condition only
            }
        }
    }

    double valueAt(int i, int j) const override {
        return field_[static_cast<size_t>(j) * static_cast<size_t>(n_) + static_cast<size_t>(i)];
    }
    std::vector<Parameter>& parameters() override { return params_; }

private:
    int n_;
    double t_ = 0.0;
    std::vector<double> field_;
    std::vector<Parameter> params_{{"kappa", 0.1, 0.0, 1.0, 0.01}};
};

} // namespace pde::examples

#pragma once

#include <string>
#include <vector>

namespace pde {

// A single tunable knob exposed as a slider in the sidebar.
struct Parameter {
    std::string name;
    double value;
    double min;
    double max;
    double step = 0.0; // 0 => continuous
};

enum class FieldKind { Scalar1D, Scalar2D };

// Implement this interface with your own numerical scheme. The GUI never
// touches the math directly -- it only calls these methods to drive the
// simulation and to read values for plotting.
class IPDEProblem {
public:
    virtual ~IPDEProblem() = default;

    virtual std::string name() const = 0;
    virtual FieldKind kind() const = 0;

    virtual int nx() const = 0;
    virtual int ny() const { return 1; }
    virtual double dx() const { return 1.0; }
    virtual double dy() const { return 1.0; }

    virtual double time() const = 0;

    // Advance the solution by dt. Put your time-marching scheme here.
    virtual void step(double dt) = 0;

    // Restore the initial condition and reset the clock.
    virtual void reset() = 0;

    // Field value at grid node (i, j). j is ignored for 1D problems.
    virtual double valueAt(int i, int j = 0) const = 0;

    // Tunable parameters (e.g. diffusivity, wave speed) rendered as sliders.
    virtual std::vector<Parameter>& parameters() = 0;

    // Optionally pin the plot's value axis/colorbar to a fixed [lo, hi]
    // range instead of auto-scaling to the current frame's min/max.
    virtual bool fixedRange(double& lo, double& hi) const {
        (void)lo;
        (void)hi;
        return false;
    }
};

} // namespace pde

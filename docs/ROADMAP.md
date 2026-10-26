# Numerical PDE — Architecture & Roadmap

Goal: solve boundary value problems for **any linear PDE up to second order**, on a
configurable domain, with configurable boundary and initial conditions.

Decisions taken (2026-08-19):

- **Linear algebra:** Eigen, pulled in via `FetchContent`. Effort goes into the
  abstract layers (operator representation, assembly, boundary conditions), not
  into reimplementing sparse factorizations.
- **Time:** steady and time-dependent problems are both first-class. Time is
  modelled explicitly rather than as "just another axis" (see §2).
- **Discretization:** finite differences on structured grids. Rectangular
  (tensor-product) domains first, then extend to the CSG shapes `Domain2d`
  already describes.

---

## 1. Where the project stands

| Area | State |
|---|---|
| `core/Expr.hpp` | **Done.** Copyable expression tree (`Number`, `Coord`, `Time`, `Binary`, `Unary`), evaluator, constant folding, operator overloads, precedence-aware printing. Replaces `expression.hpp` / `.cpp`. |
| `core/MultiIndex.hpp` | **Done.** Derivative multi-index. Mixed partials canonicalise themselves; terms sort principal-part-first. |
| `core/Operator.hpp` | **Done.** `Term`, `LinearOperator`, `PDE` in the canonical form of §2, plus builders for Laplace / Poisson / Helmholtz / heat / wave / advection–diffusion. Replaces `PDE_term.hpp`, `GeneralizedFirstOrder.hpp`, `GeneralizedSecondOrder.hpp`. |
| `Domain1d.hpp` | Works for a single interval. |
| `Domain2d.hpp/.cpp` | CSG add/subtract of rectangles and circles. Several correctness bugs; `inDomain` has no `return`. |
| `meshes/mesh1d.hpp`, `mesh2d.hpp` | Empty stubs, neither compiles. |
| `gui/` | ~700 lines of FLTK app (`AppWindow`, `PlotWidget`, `ColorMap`) behind a clean `IPDEProblem` interface. **The most finished part of the project.** Its `CMakeLists.txt` is stale — see appendix. |
| Build system | None yet. `core/` is header-only and compiles standalone. |
| Tests | `tests/test_operator.cpp`, run by hand until Phase 0 wires up CTest. |

The GUI already draws the right boundary: `IPDEProblem` exposes `step(dt)`,
`reset()`, `valueAt(i,j)`, and tunable parameters, and knows nothing about the
math. That contract is worth keeping exactly as it is — the solver becomes an
implementation of it.

The original structural problem — **three competing designs** for the same thing
(`GeneralizedFirstOrder`, named coefficients plus time; `GeneralizedSecondOrder`,
named coefficients without time; `PDE_Term`, coefficient plus order plus
derivative factors) — is **resolved**. Only the third generalized, and §2/§3
built it out; the other two are deleted. `tests/test_operator.cpp` includes
worked mappings showing both old structs expressed in the new form.

---

## 2. One canonical form

Every linear PDE up to second order, in any number of spatial dimensions, fits:

$$
m_2(\mathbf{x})\,\partial_t^2 u \;+\; m_1(\mathbf{x})\,\partial_t u
\;=\;
\sum_{i \le j} a_{ij}(\mathbf{x},t)\,\partial_i \partial_j u
\;+\; \sum_i b_i(\mathbf{x},t)\,\partial_i u
\;+\; c(\mathbf{x},t)\,u
\;+\; f(\mathbf{x},t)
$$

Spatial operator on the right, time derivatives on the left. This single form
covers the whole target space:

| Problem | Coefficients |
|---|---|
| Laplace / Poisson | `m₂=m₁=0`, `a_ii=1`, `f=−ρ` |
| Helmholtz | `m₂=m₁=0`, `a_ii=1`, `c=k²` |
| Heat / diffusion | `m₁=1`, `a_ii=κ` |
| Advection–diffusion | `m₁=1`, `a_ii=κ`, `b_i=−v_i` |
| Wave | `m₂=1`, `a_ii=c²` |
| Anisotropic / variable-coefficient anything | `a_ij`, `b_i`, `c` as arbitrary expressions of `x` and `t` |

`GeneralizedFirstOrder` and `GeneralizedSecondOrder` both collapse into this.
**Delete both** once §3 lands.

Why time is separated rather than treated as a fourth axis: a BVP is elliptic in
space but an *initial* value problem in time. Making `t` an axis would demand
boundary data at the final time, which nobody has. Keeping `m₂`/`m₁` on the left
means one assembled spatial operator serves both the steady solve and every time
step.

---

## 3. The representation that makes this one code path

Do not name coefficients. Key each term by a **derivative multi-index**.

```cpp
// include/pde/core/MultiIndex.hpp
// alpha[i] = order of differentiation along axis i.  sum(alpha) <= 2.
//   u          -> {0,0}
//   du/dx      -> {1,0}
//   d2u/dxdy   -> {1,1}
//   d2u/dy2    -> {0,2}
using MultiIndex = std::array<uint8_t, kMaxDim>;

// include/pde/core/Operator.hpp
struct Term {
    MultiIndex alpha;
    Expr       coefficient;   // function of x_0..x_{d-1}, t
};

class LinearOperator {          // the right-hand side of §2
    int dim_;
    std::vector<Term> terms_;
    Expr source_;               // f
public:
    LinearOperator& add(MultiIndex alpha, Expr coeff);
    // convenience builders
    static LinearOperator laplacian(int dim);
    static LinearOperator advectionDiffusion(int dim, Expr kappa, std::vector<Expr> v);
};

struct PDE {
    LinearOperator spatial;
    Expr mass1{Number{0.0}};    // coefficient of du/dt
    Expr mass2{Number{0.0}};    // coefficient of d2u/dt2
};
```

This is where `PDE_term.hpp` was already heading — replace its
`std::pair<std::string,std::string> factor1/factor2` with a `MultiIndex`. The
multi-index is comparable, hashable, dimension-agnostic, and maps directly onto
a stencil, which the string pairs never could.

### The key insight: stencils are tensor products

You only ever need **three one-dimensional stencils**:

```
order 0 (identity):   { (0, 1) }                        weights at offsets
order 1 (central):    { (-1, -1/2h), (+1, +1/2h) }
order 2 (central):    { (-1, 1/h²), (0, -2/h²), (+1, 1/h²) }
```

The stencil for **any** multi-index is their tensor product across axes:

$$ S(\alpha) = \bigotimes_{i=0}^{d-1} S_{1D}(\alpha_i, h_i) $$

Worked example — the mixed derivative `∂²u/∂x∂y` falls out automatically as
`S_1D(1,hx) ⊗ S_1D(1,hy)`:

```
(+1,+1): +1/(4 hx hy)      (+1,-1): -1/(4 hx hy)
(-1,+1): -1/(4 hx hy)      (-1,-1): +1/(4 hx hy)
```

which is the standard second-order mixed stencil. Consequences:

- 2D, 3D, and nD need **no new code** — only a longer multi-index.
- Higher-order accuracy later = swapping the 1D stencil table, nothing else.
- Fourth-order PDEs later = adding `order 3` and `order 4` rows to that table.

This is the abstraction that makes "any PDE up to second order" genuinely a
single code path rather than a switch over equation types. Build it early.

---

## 4. Module layout

```
CMakeLists.txt                          top-level; Eigen + FLTK via FetchContent
include/pde/
  core/Expr.hpp                         AST, eval, symbolic diff, (later) parser
  core/MultiIndex.hpp
  core/Operator.hpp                     Term, LinearOperator, PDE
  geometry/Region.hpp                   d-dimensional shapes + CSG
  geometry/Domain.hpp                   inside / on-boundary / outward normal
  mesh/StructuredGrid.hpp               spacing, node classification, DOF map
  discretize/Stencil.hpp                1D table + tensor product
  discretize/BoundaryCondition.hpp      Robin form, per-segment assignment
  assemble/Assembler.hpp                PDE + Grid + BCs -> sparse system
  solve/LinearSolver.hpp                Eigen wrapper behind an interface
  solve/SteadySolver.hpp
  solve/TimeStepper.hpp                 theta-method, first-order reduction
  verify/Manufactured.hpp               MMS harness
src/                                    .cpp for non-template code
tests/
apps/gui/                               existing FLTK app, moved from src/gui
examples/
```

Dependency direction is strictly downward: `core` knows nothing about meshes,
`mesh` knows nothing about assembly, `gui` depends on everything and is depended
on by nothing.

---

## 5. Phases

Each phase ends with something runnable. Do not start a phase before the
previous one's "done when" holds.

### Phase 0 — Unblock (half a day)

Nothing below can proceed while the tree does not compile.

1. Top-level `CMakeLists.txt`: target `pde_core` (header-only INTERFACE or a
   small static lib), `FetchContent` Eigen, `enable_testing()`.
2. Fix the seven compile errors and the duplicate include guard (appendix).
3. Delete `expression.cpp`; keep `expression.hpp` as the single definition.
4. Clean up the git index — the staged `generalized1d.hpp -> generalized_first_order.hpp`
   renames are half-applied and the CamelCase replacements are untracked.
5. Add a test target with one trivial passing test, so CI has a shape.

**Done when:** `cmake --build build && ctest` is green from a clean clone.

### Phase 1 — Expression layer (mostly done)

`Expr` is load-bearing: it carries coefficients, sources, boundary data, initial
conditions, and (via §7) manufactured solutions.

- [x] Ops beyond add/multiply. Collapsed into one `Binary{BinOp}` node plus one
      `Unary{UnFn}` node rather than a struct per operation — adding an
      operation is now one enum value and one switch case, not a new variant
      alternative and a new branch in every traversal (`eval`, `toString`, and
      later `diff` and `simplify` are each a full traversal).
- [x] `eval(const Expr&, const Env&)` in the library, with `Env` holding
      `{const double* x; int dim; double t;}`. Coordinates are addressed by
      index (`Coord{axis}`, `Time{}`) instead of by name, so evaluation is an
      array load rather than a string comparison — it runs once per matrix entry
      during assembly. Out-of-range axes throw rather than returning a silent
      zero.
- [x] Value semantics, via `shared_ptr` on the interior nodes. Trees are
      immutable, so sharing is safe and a copy is a pointer bump. This is what
      lets a whole `PDE` be a copyable value.
- [x] Constant folding in the constructors, plus precedence-aware printing.
- [ ] **Symbolic differentiation** — `Expr diffAxis(const Expr&, int axis)` and
      `Expr diffTime(const Expr&)`. Now the only substantial piece left, and the
      one that makes Phase 7 free. Note the signature follows the indexed
      coordinate scheme rather than taking a variable name.

**Done when:** tests cover `diff` against hand-computed derivatives.

*Optional, deferrable:* a string parser (`"sin(x)*exp(-t)"`). Nice for the GUI
and config files, but it is a self-contained side quest — do not let it block
the solver.

### Phase 2 — Operator representation ✅ done

§3 implemented in `core/MultiIndex.hpp` and `core/Operator.hpp`:
`MultiIndex`, `Term`, `LinearOperator`, `PDE`, `TimeKind`, and builders for
Laplace / Poisson / Helmholtz / heat / wave / advection–diffusion.
`GeneralizedFirstOrder.hpp`, `GeneralizedSecondOrder.hpp`, and `PDE_term.hpp`
are deleted.

Design points worth remembering:

- `add()` **accumulates** on a repeated multi-index rather than duplicating, and
  drops the term when coefficients cancel. The term list is therefore canonical:
  at most one entry per derivative, so assembly never walks a zero term.
- Terms are stored principal-part-first (`u_xx`, `u_xy`, `u_yy`, `u_x`, `u_y`,
  `u`), which makes the printed operator read like the textbook form.
- `TimeKind` classifies the **time structure** only (`Steady` / `FirstOrder` /
  `SecondOrder`) because that is what selects the solver path. It deliberately
  does not claim ellipticity/parabolicity/hyperbolicity — that needs the
  eigenvalues of the principal part and is pointwise when coefficients vary.
- There is intentionally no `general2d(A, B, C, D, E, F)` builder; that would
  reintroduce the naming problem this form removes.

**Done:** heat and wave are one line each, `tests/test_operator.cpp` asserts
their term lists, and it carries worked mappings of both retired structs into
the new form.

### Phase 3 — Grid and degrees of freedom (1–2 days)

```cpp
class StructuredGrid {
    int dim_;
    std::array<int, kMaxDim>    n_;      // nodes per axis
    std::array<double, kMaxDim> origin_, h_;
public:
    enum class NodeType { Interior, Boundary, Exterior };
    // flat index <-> multi-index
    int  flatten(const std::array<int,kMaxDim>&) const;
    void unflatten(int, std::array<int,kMaxDim>&) const;
    void coords(int flat, double* x) const;

    NodeType classify(int flat) const;
    int      dof(int flat) const;        // -1 for Exterior
    int      numDofs() const;
};
```

Boundary nodes **keep** their DOF — their matrix row becomes the boundary
condition. Eliminating them looks tidier but makes Neumann and Robin awkward and
makes the residual harder to inspect.

Node classification is the seam where arbitrary geometry plugs in later
(Phase 9): for now, `classify` reads the rectangle; later it queries `Domain`.

**Done when:** classification and DOF numbering are tested on 1D, 2D, and 3D
grids, including a grid with an interior hole punched out (exterior nodes
correctly skipped even though nothing solves on it yet).

### Phase 4 — Stencils (1 day)

The tensor-product machinery from §3.

```cpp
struct StencilEntry { std::array<int,kMaxDim> offset; double weight; };
using Stencil = std::vector<StencilEntry>;

Stencil stencilFor(MultiIndex alpha, const double* h, int dim);
```

**Done when:** tests verify `∂²/∂x²`, `∂/∂x`, and `∂²/∂x∂y` against the known
weights above, **and** a stencil applied to a polynomial reproduces its exact
derivative (the sharpest test available — central differences are exact on
low-degree polynomials).

### Phase 5 — Boundary conditions and assembly (2–4 days)

The heart of the project.

One boundary condition type covers all three classical kinds — the Robin form:

$$ \alpha(\mathbf{x})\,u + \beta(\mathbf{x})\,\partial_n u = g(\mathbf{x},t) $$

`(1, 0, g)` is Dirichlet, `(0, 1, g)` is Neumann, general is Robin. Assign them
by predicate or by named boundary segment, with a default fallback.

```cpp
struct BoundaryCondition { Expr alpha, beta, g; };

class Assembler {
public:
    struct System { Eigen::SparseMatrix<double> A; Eigen::VectorXd b; };
    System assembleSteady(const PDE&, const StructuredGrid&, const BCSet&, double t = 0.0);
    Eigen::SparseMatrix<double> assembleMass(const Expr& m, const StructuredGrid&);
};
```

Assembly loop, for each node with a DOF:

- **Interior:** for each `Term`, evaluate its coefficient at the node, fetch
  `stencilFor(alpha, h, dim)`, scatter `coeff * weight` into the row. Source `f`
  goes to `b`.
- **Boundary:** emit the BC row instead. `∂_n u` on an axis-aligned face uses a
  second-order one-sided difference, e.g. at a left-`x` face
  `∂_n u = -(-3u_0 + 4u_1 - u_2)/(2h)`.

Build with `std::vector<Eigen::Triplet<double>>` then `setFromTriplets` — do not
insert into a `SparseMatrix` element by element.

Watch for: a stencil reaching an exterior node (must be caught, not silently
dropped); the sign convention on outward normals; and rows that end up
structurally empty.

**Done when:** assembling 1D Laplace with Dirichlet ends produces exactly the
tridiagonal `(1,−2,1)/h²` matrix, verified entry by entry against a
hand-written expected matrix.

### Phase 6 — Steady solve (half a day)

Thin wrapper. Keep it behind an interface so the solver is swappable:

```cpp
class LinearSolver { public: virtual Eigen::VectorXd solve(const SparseMatrix&, const VectorXd&) = 0; };
```

Default to `Eigen::SparseLU` (robust, direct, fine to ~10⁵ unknowns). Add
`BiCGSTAB` + `IncompleteLUT` for large or 3D problems. `ConjugateGradient` only
for the symmetric case — note that Dirichlet rows and any first-derivative term
destroy symmetry, so do not reach for CG by default.

**Done when:** 1D Poisson with `u(0)=u(1)=0`, `f=−1` returns `x(1−x)/2` to
round-off, and 2D Laplace on a square with a nonconstant Dirichlet boundary
matches a separable analytic solution.

### Phase 7 — Verification harness (1 day)

This is the payoff for symbolic differentiation, and the thing that will tell
you whether everything above is actually right.

**Method of manufactured solutions.** Choose any `u_exact` as an `Expr`. Apply
the operator to it *symbolically* using `diff` to obtain the exact `f` and the
exact boundary data. Solve. Compare. Refine the grid and check the error ratio.

```cpp
Expr applySymbolically(const LinearOperator&, const Expr& u);   // uses diff()
double convergenceOrder(const PDE&, const Expr& u_exact, std::vector<int> resolutions);
```

Second-order central differences must give an observed order of ≈2.0. An order
of 1.0 means a boundary condition is only first-order accurate; an order of 0
means something is plainly wrong. **This single test catches more bugs than
every other test combined.** Wire it into `ctest` over a spread of manufactured
solutions and coefficient sets.

**Done when:** observed order is 2.0 ± 0.1 for Poisson, Helmholtz, anisotropic
diffusion (nonzero `a_xy`), and advection–diffusion, in both 1D and 2D.

### Phase 8 — Time stepping (2–3 days)

The spatial operator `L` and mass matrices are already assembled. If
coefficients are time-independent, assemble **once** and reuse across all steps
— and factorize once too, for implicit schemes with fixed `Δt`. That single
optimization is the difference between an interactive GUI and a slideshow.

**Parabolic** (`m₂ = 0`) — θ-method:

$$ (M_1 - \theta \Delta t\, L)\,u^{n+1} = (M_1 + (1-\theta)\Delta t\, L)\,u^n + \Delta t\left[\theta f^{n+1} + (1-\theta) f^n\right] $$

`θ=0` explicit Euler, `θ=½` Crank–Nicolson (second order, the default),
`θ=1` backward Euler (first order, maximally stable).

**Hyperbolic** (`m₂ ≠ 0`) — reduce to a first-order system with `v = ∂u/∂t`:

$$ \begin{bmatrix} I & 0 \\ 0 & M_2 \end{bmatrix} \frac{d}{dt}\begin{bmatrix} u \\ v \end{bmatrix} = \begin{bmatrix} 0 & I \\ L & -M_1 \end{bmatrix}\begin{bmatrix} u \\ v \end{bmatrix} + \begin{bmatrix} 0 \\ f \end{bmatrix} $$

then apply the *same* θ-method to the block system. One integrator, both
equation classes. (Do not route the `m₂=0` case through the block form — it
degenerates into a DAE. Branch on `m₂ ≡ 0`.)

Explicit schemes need a stability guard, computed and reported rather than
discovered by watching the solution explode:

- diffusion: `Δt ≤ h² / (2 d κ)`
- advection: `Δt ≤ h / |v|`
- wave: `Δt ≤ h / (c √d)`

Initial conditions: `u(x,0)` as an `Expr`; a second `Expr` for `∂u/∂t(x,0)` when
`m₂ ≠ 0`.

**Done when:** heat equation on `[0,1]` with a sine initial condition decays at
the analytic rate `exp(−κπ²t)`; the wave equation conserves energy to within
scheme tolerance over many periods; MMS with a time-dependent exact solution
shows second order in *both* `h` and `Δt` for Crank–Nicolson.

### Phase 9 — Arbitrary geometry (3–5 days)

Only now does the `Domain2d` CSG work pay off.

1. Fix and generalize the geometry (appendix has the bug list). Make it
   d-dimensional: `Region` with `Box`, `Ball`, `HalfSpace`, and a CSG tree of
   union / intersection / difference. Add `signedDistance` and
   `outwardNormal` — both are needed for cut-cell boundary treatment and neither
   is expressible with the current `contains`/`on_boundary` pair alone.
2. Point `StructuredGrid::classify` at the `Domain`.
3. Handle boundaries that fall *between* grid nodes — Shortley–Weller: when the
   boundary cuts the link between nodes at fraction `θ ∈ (0,1)`, the stencil
   weights adjust to the uneven spacing. Start with Dirichlet only; Neumann on a
   cut cell is materially harder.

**Done when:** Laplace on an annulus (square minus a disc) with different
Dirichlet values on the inner and outer boundary matches the analytic
logarithmic solution, and MMS on a disc still shows order ≥ 1.5.

### Phase 10 — Wire up the GUI (1 day)

Write one adapter and the existing app comes alive unchanged:

```cpp
class SolverProblem : public IPDEProblem {
    // holds PDE, StructuredGrid, BCSet, TimeStepper, current solution
    void   step(double dt) override;        // one time step, or re-solve steady
    void   reset() override;                // re-apply initial condition
    double valueAt(int i, int j) const override;
    std::vector<Parameter>& parameters() override;   // bind to Expr constants
};
```

Interesting detail worth designing for: `parameters()` exposes sliders, so
binding a named constant in an `Expr` to a slider gives live re-solve as the
user drags `κ`. That requires assembly to be fast enough — hence the
assemble-once caching in Phase 8, and cache invalidation when a bound parameter
changes.

Also replace `TemplateProblems.hpp` with real examples, and fix the stale
`CMakeLists.txt` (appendix).

**Done when:** the FLTK app animates a real heat-equation solution and dragging
`κ` visibly changes the decay rate.

---

## 6. Beyond

Roughly in order of value:

- **Nonlinear terms** — coefficients depending on `u`. Needs Newton iteration
  and a Jacobian; the assembler already computes almost the right thing.
- **Systems of PDEs** — vector `u`. Multiply DOF count by the number of
  components and generalize `Term` to carry `(row, col)` component indices.
- **Higher-order stencils** — one new table, per §3.
- **Adaptive refinement** — a real change; structured grids fight back. Consider
  quadtree/octree only if genuinely needed.
- **Parser + config files** — define problems in text instead of C++.
- **Parallel assembly** — the loop is embarrassingly parallel per row.

Two things worth resisting for now: templating on dimension
(`Grid<2>`, `Grid<3>`) — runtime `dim_` with a `kMaxDim` cap is far simpler and
the cost is negligible; and abstracting the discretization behind an interface
so FEM could slot in later — that interface cannot be designed correctly before
having written one FEM assembler, and a premature version will be wrong.

---

## Appendix — bug inventory (Phase 0 checklist)

Confirmed by compiling each header with `g++ -std=c++17 -fsyntax-only`.

**Resolved by Phase 1/2** — the three files that did not compile
(`PDE_term.hpp`, `GeneralizedFirstOrder.hpp`, `GeneralizedSecondOrder.hpp`) are
deleted, along with the duplicate `expression.hpp` / `expression.cpp` pair that
both defined `pde::Expr` in the same namespace. That also removes the
`"Expression.hpp"`-with-a-capital-E includes, which worked only because Windows
filesystems are case-insensitive. `src/core/` compiles clean under
`-Wall -Wextra -Wpedantic`.

Note that `_scratch/test2.cpp` included `src/expression.cpp` directly and will
no longer build; its `eval` demo now lives in the library and is covered by
`tests/test_operator.cpp`.

**Does not compile**

- `src/meshes/mesh1d.hpp` — `#include "Domain1d.hpp"` should be
  `"../Domain1d.hpp"`; missing `;` after the `struct` closing brace.
- `src/meshes/mesh2d.hpp` — includes `"Domain.hpp"`, which no longer exists
  (deleted in `fc49ded`); missing `;` after the `struct`.

**Silent, worse than a compile error**

- `Domain1d.hpp` and `Domain2d.hpp` **both use `#ifndef DOMAIN_HPP`**. Including
  both in one translation unit silently drops the second — verified: after
  including `Domain1d.hpp`, `pde::Domain2d` does not exist. Give each a unique
  guard, or use `#pragma once` as the GUI headers already do.
- Both also define `constexpr double EPSILON` in namespace `pde`; once the
  guards are fixed, that becomes a redefinition. Hoist it to one shared header
  — `src/core/MultiIndex.hpp` already holds `kMaxDim` and is the natural home.
- `Domain2d::inDomain` (`Domain2d.cpp:7`) computes `is_inside` and then **falls
  off the end without returning it** — undefined behaviour.

**Wrong math**

- `Circle::contains` (`Domain2d.hpp:22`) compares squared distance against
  `radius`, not `radius * radius`.
- `Circle::on_boundary` uses exact `==` on doubles; effectively never true. Use
  the `EPSILON` tolerance that `Rectangle::on_boundary` already uses.
- `Domain2d::inDomain` subtraction: `is_inside = is_inside && inOperation`
  should be `is_inside && !inOperation`. As written, subtract behaves as
  intersect.
- `Domain2d::onBoundary` computes `onOperationBoundary` and never reads it; the
  branches only ever assign `false` or leave the value unchanged, so the
  function always returns `false`. Needs rewriting — a point is on the CSG
  boundary if it lies on some primitive's boundary *and* that primitive's
  surface survives the remaining operations. This is the single subtlest piece
  of the geometry layer; a signed-distance formulation (`min`/`max`/negate) is
  much easier to get right than boolean bookkeeping.
- `Domain1d::inDomain(double x, double y)` takes an unused `y` — signature
  copy-pasted from 2D.

**Build**

- `src/gui/CMakeLists.txt` declares `add_executable(pde_solver main.cpp)`, but
  the committed build artifacts are for a target `pde_studio` built from
  `src/main.cpp`, `src/AppWindow.cpp`, `src/PlotWidget.cpp`. The file is stale
  and **cannot configure from a clean checkout** (no `main.cpp` at that path, no
  include directories, wrong target name). It also `FetchContent`s FLTK from
  git, while the built binary links MSYS2's system FLTK
  (`C:/msys64/ucrt64/lib/libfltk*.dll.a`). Pick one and make the file match.
- No build system exists for the core library at all — `.vscode/tasks.json`
  compiles single files with a hardcoded `g++.exe` path.

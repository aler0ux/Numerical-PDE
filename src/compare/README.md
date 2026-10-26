# FEM vs FDM comparison (Poisson equation)

A small, self-contained learning tool: solve

    -Laplacian(u) = f     on a chosen domain,     u = 0 on the boundary

twice -- once with a finite-difference 5-point stencil, once with linear
(P1) finite elements -- on the *same* structured grid and degree-of-freedom
numbering, so the two solutions can be compared node by node. Both
assemblies produce a plain CSR matrix and are solved with the same
hand-written Conjugate Gradient iterative solver.

## Files

- `CSRMatrix.hpp` -- Compressed Sparse Row matrix built from triplets, with `matvec`.
- `ConjugateGradient.hpp` -- Jacobi-preconditioned CG, the one iterative solver used for both methods.
- `Domains.hpp` -- four preset domains: `square`, `lshape`, `plus`, `circle`.
- `SourceFunctions.hpp` -- a small library of right-hand sides `f(x,y)`, written as `pde::Expr` (the project's existing symbolic expression type from `core/Expr.hpp`), so picking `f` reuses the project's own PDE framework instead of a separate mechanism.
- `Grid2d.hpp` -- the shared (N+1)x(N+1) structured grid, classified against a domain and numbered into DOFs.
- `Assemble.hpp` -- the two assemblers (`assembleFDM`, `assembleFEM`) plus the shared Dirichlet-elimination pass that keeps both systems symmetric positive definite.
- `main.cpp` -- CLI / interactive driver: builds both systems, solves both, times and reports both, writes a CSV.
- `plot_compare.py` -- reads the CSV and draws a 2x2 comparison figure (FDM field, FEM field, |difference|, CG convergence).
- `build.sh`, `run_compare.sh` -- build and build+solve+plot in one step.

## Usage

```sh
cd src/compare
./run_compare.sh circle sine 64
./run_compare.sh lshape bump 48
./run_compare.sh plus constant 80
./run_compare.sh square checker 100
```

Domains: `square`, `lshape`, `plus`, `circle`.
Source functions: `constant`, `sine`, `bump`, `linear`, `checker` (see
`SourceFunctions.hpp` for the exact formulas -- it is easy to add more).

Running `compare2d` with no arguments drops into an interactive menu that
prompts for the same four choices.

Each run prints, for both methods: number of unknowns, nonzero count,
CG iteration count, final residual, and wall-clock solve time, then the
maximum pointwise difference between the two solutions. It writes
`<output>.csv` and `<output>.csv.residuals.csv`, and `plot_compare.py`
turns those into `<output>.png`.

## The numerics, briefly

- **Shared grid.** Both methods run on identical nodes, so "compare FEM and
  FDM" means literally differencing the same vector of node values, with no
  interpolation between meshes.
- **Domain classification.** A node is *active* if it lies in the domain
  (`Domains.hpp::inDomain`), and *boundary* if it is active but has a
  grid-edge or inactive neighbour -- i.e. where the discrete domain meets
  its own edge. This is the simplest ("immersed"/staircase) way to put an
  arbitrary shape on a regular grid: the boundary is only first-order
  accurate on a curved edge (the circle), which is the expected tradeoff
  for keeping the geometry this simple.
- **FDM.** The standard 5-point stencil, `(4u_0 - sum(neighbours)) / h^2 =
  f`, assembled only at non-boundary active nodes (which are guaranteed to
  have all four neighbours active).
- **FEM.** Each grid cell with >= 3 active corners becomes one or two
  linear triangles (the usual structured-mesh diagonal split when all four
  corners are active). Element stiffness uses the closed-form P1 formula
  `K_ij = (b_i b_j + c_i c_j) / (4*Area)`; the load vector uses one-point
  (centroid) quadrature, exact for constant `f` and a standard simple
  approximation otherwise.
- **Boundary condition.** Homogeneous Dirichlet (`u = 0`) is applied by
  zeroing every row *and column* of a boundary DOF and placing a `1` on the
  diagonal (`Assemble.hpp::applyDirichletZero`). Because the boundary value
  is zero there is no right-hand-side correction to propagate -- dropping
  those columns is exactly that correction -- and the result stays
  symmetric positive definite, so Conjugate Gradient applies directly to
  both systems without special-casing either one.
- **Solver.** Jacobi-preconditioned CG on the CSR matrix, convergence
  tracked as the relative residual `||b - Ax|| / ||b||`.

## Deliberately out of scope

This lives beside `docs/ROADMAP.md`'s longer-term plan (Eigen-backed
assembler, general CSG geometry via `Domain2d`, time-dependent problems,
the GUI) rather than inside it: it is a focused, from-scratch, CSR +
hand-rolled-CG exercise for comparing two discretizations on a few fixed
domains, for one equation (Poisson). It does not touch `Domain2d.hpp`,
`mesh2d.hpp`, or the GUI, and does not depend on Eigen.

## Interactive desktop window (`src/compare/gui`)

A standalone FLTK window for exploring the comparison live, instead of
through the CLI + Python plot: pick the domain, `f`, the mesh resolution
`N`, and the method (FDM only / FEM only / both side by side) from the
sidebar, then press **Solve** to re-assemble and re-run Conjugate Gradient
and refresh the field(s) on the right. Nothing re-solves automatically
when a control changes -- Solve is the explicit "refresh" step.

A separate **View** dropdown switches between a flat **Heatmap** and a
shaded **3D surface** (a software axonometric/isometric projection, no
OpenGL dependency) for whatever was last solved -- this one *does* switch
immediately, since it's just a different rendering of the same numbers,
not a re-solve.

Files (`src/compare/gui/`):
- `CompareWindow.hpp/.cpp` -- the window: sidebar controls + the solve-and-refresh logic, calling the exact same `assembleFDM`/`assembleFEM`/`conjugateGradient` used by the CLI.
- `HeatmapView.hpp/.cpp` -- a small `Fl_Widget` that draws one or two node-value heatmaps with colorbars, reusing `pde::viridis` from `src/gui/include/ColorMap.hpp` (the same colormap the existing `pde_studio` app uses).
- `main.cpp`, `CMakeLists.txt` -- entry point and build file for this target (`compare_gui`), independent of `src/gui`'s own (stale) `CMakeLists.txt`.

### Building

```sh
cd src/compare/gui
cmake -S . -B build -G "MinGW Makefiles"   # or your usual generator
cmake --build build
./build/compare_gui
```

The `CMakeLists.txt` first looks for an installed FLTK (`find_package(FLTK)`
-- on MSYS2, `pacman -S mingw-w64-ucrt-x86_64-fltk`) and only falls back to
fetching and building FLTK from source via `FetchContent` if none is found,
since the latter is slow on a first build.

> This GUI was written and carefully cross-checked against the existing,
> already-building `src/gui/src/AppWindow.cpp` / `PlotWidget.cpp` (same FLTK
> widget and drawing idioms throughout), but it has **not been compiled** in
> this environment -- the sandbox used to build and test the rest of this
> module has no FLTK installation and no permission to build one from
> source here. Build it on your machine and let me know what, if anything,
> doesn't compile.

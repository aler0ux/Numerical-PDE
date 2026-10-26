// FEM vs FDM comparison driver for the Poisson equation
//     -Laplacian(u) = f   on a chosen domain,   u = 0 on the boundary
//
// Usage:
//     compare2d [domain] [source] [N] [output.csv]
//     compare2d                          interactive menu
//
//   domain  : square | lshape | plus | circle | eight | donut  (default: circle)
//   source  : constant | sine | bump | linear | checker (default: sine)
//   N       : grid cells per axis, grid is (N+1)x(N+1) nodes (default: 64)
//   output  : CSV path the Python plotting script reads (default: result.csv)
//
// Solves the same problem twice -- once with finite differences, once with
// linear finite elements -- on the same Grid2d, using the same CSR +
// Conjugate-Gradient machinery for both, and writes a single CSV with both
// node fields side by side so plot_compare.py can render them.

#include <algorithm>
#include <chrono>
#include <cmath>
#include <fstream>
#include <iostream>
#include <string>

#include "Assemble.hpp"
#include "ConjugateGradient.hpp"
#include "Domains.hpp"
#include "Grid2d.hpp"
#include "SourceFunctions.hpp"
#include "../core/Operator.hpp"

namespace {

std::string prompt(const std::string& label, const std::string& def) {
    std::cout << label << " [" << def << "]: ";
    std::string line;
    std::getline(std::cin, line);
    if (line.empty()) return def;
    return line;
}

void printMenu() {
    std::cout << "Domains: ";
    for (auto k : {cmp::DomainKind::Square, cmp::DomainKind::LShape, cmp::DomainKind::Plus, cmp::DomainKind::Circle,
                   cmp::DomainKind::Eight, cmp::DomainKind::Donut})
        std::cout << cmp::domainName(k) << "  ";
    std::cout << "\nSource functions:\n";
    for (const auto& f : cmp::sourceFunctionLibrary()) std::cout << "  " << f.name << " -- " << f.description << "\n";
    std::cout << "\n";
}

}  // namespace

int main(int argc, char** argv) {
    std::string domainArg = "circle";
    std::string sourceArg = "sine";
    int n = 64;
    std::string outputPath = "result.csv";

    if (argc >= 2) {
        domainArg = argv[1];
        sourceArg = (argc >= 3) ? argv[2] : sourceArg;
        n = (argc >= 4) ? std::stoi(argv[3]) : n;
        outputPath = (argc >= 5) ? argv[4] : outputPath;
    } else {
        printMenu();
        domainArg = prompt("Domain", domainArg);
        sourceArg = prompt("Source function", sourceArg);
        n = std::stoi(prompt("Grid cells per axis N", std::to_string(n)));
        outputPath = prompt("Output CSV path", outputPath);
    }

    const cmp::DomainKind domain = cmp::domainFromName(domainArg);
    const cmp::SourceFunction& source = cmp::sourceFunctionFromName(sourceArg);

    pde::PDE problem = pde::PDE::poisson(2, source.expr);
    std::cout << "\nSolving:  " << problem.toString() << "\n";
    std::cout << "Domain:   " << cmp::domainName(domain) << "\n";
    std::cout << "Grid:     " << (n + 1) << " x " << (n + 1) << " nodes (N=" << n << ")\n\n";

    cmp::Grid2d grid(domain, n);
    std::cout << "Active nodes (shared DOF count): " << grid.numDofs() << "\n\n";

    auto solveAndReport = [](const char* label, const cmp::System& sys) {
        const auto t0 = std::chrono::steady_clock::now();
        cmp::CGResult result = cmp::conjugateGradient(sys.A, sys.b);
        const auto t1 = std::chrono::steady_clock::now();
        const double ms = std::chrono::duration<double, std::milli>(t1 - t0).count();
        std::cout << label << ":\n"
                  << "  unknowns     = " << sys.A.rows() << "\n"
                  << "  nonzeros     = " << sys.A.nnz() << "\n"
                  << "  CG iters     = " << result.iterations << (result.converged ? " (converged)" : " (NOT converged)") << "\n"
                  << "  final resid  = " << result.residualHistory.back() << "\n"
                  << "  solve time   = " << ms << " ms\n\n";
        return result;
    };

    const cmp::System fdmSystem = cmp::assembleFDM(grid, source.expr);
    cmp::CGResult fdmResult = solveAndReport("Finite differences (5-point stencil)", fdmSystem);

    const cmp::System femSystem = cmp::assembleFEM(grid, source.expr);
    cmp::CGResult femResult = solveAndReport("Finite elements (P1 triangles)", femSystem);

    double maxDiff = 0.0;
    for (int d = 0; d < grid.numDofs(); ++d) maxDiff = std::max(maxDiff, std::fabs(fdmResult.x[static_cast<size_t>(d)] - femResult.x[static_cast<size_t>(d)]));
    std::cout << "Max |u_fdm - u_fem| over shared nodes: " << maxDiff << "\n\n";

    // --- write the CSV the plotting script consumes -------------------
    std::ofstream out(outputPath);
    if (!out) {
        std::cerr << "Could not open output file '" << outputPath << "' for writing\n";
        return 1;
    }
    out << "# domain=" << cmp::domainName(domain) << " source=" << source.name << " n=" << n
        << " fdm_iters=" << fdmResult.iterations << " fem_iters=" << femResult.iterations << "\n";
    out << "i,j,x,y,active,boundary,u_fdm,u_fem\n";
    for (int i = 0; i <= grid.n(); ++i) {
        for (int j = 0; j <= grid.n(); ++j) {
            const bool act = grid.active(i, j);
            const double ufdm = act ? fdmResult.x[static_cast<size_t>(grid.dof(i, j))] : 0.0;
            const double ufem = act ? femResult.x[static_cast<size_t>(grid.dof(i, j))] : 0.0;
            out << i << "," << j << "," << grid.x(i) << "," << grid.y(j) << "," << (act ? 1 : 0) << ","
                << (act && grid.boundary(i, j) ? 1 : 0) << "," << ufdm << "," << ufem << "\n";
        }
    }

    // A second small file for the residual-history plot (CG iteration curves).
    std::ofstream resFile(outputPath + ".residuals.csv");
    resFile << "iteration,fdm_residual,fem_residual\n";
    const size_t maxIters = std::max(fdmResult.residualHistory.size(), femResult.residualHistory.size());
    for (size_t k = 0; k < maxIters; ++k) {
        const double rf = (k < fdmResult.residualHistory.size()) ? fdmResult.residualHistory[k] : fdmResult.residualHistory.back();
        const double rm = (k < femResult.residualHistory.size()) ? femResult.residualHistory[k] : femResult.residualHistory.back();
        resFile << (k + 1) << "," << rf << "," << rm << "\n";
    }

    std::cout << "Wrote " << outputPath << " and " << outputPath << ".residuals.csv\n";
    return 0;
}

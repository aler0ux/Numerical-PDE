#!/usr/bin/env bash
# Build the FEM-vs-FDM comparison executable.
#
# On the project's normal MSYS2/UCRT64 toolchain, run this from an MSYS2
# shell (or adapt the compiler path) -- it needs nothing but a C++17
# compiler, no external libraries.
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")"
g++ -std=c++17 -O2 -Wall -Wextra -o compare2d main.cpp
echo "Built ./compare2d"

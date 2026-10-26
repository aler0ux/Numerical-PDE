#!/usr/bin/env bash
# Build (if needed), solve, and plot in one step:
#   ./run_compare.sh [domain] [source] [N]
# e.g.  ./run_compare.sh lshape bump 48
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")"
if [ ! -x ./compare2d ] || [ main.cpp -nt ./compare2d ]; then
    ./build.sh
fi
domain="${1:-circle}"
source_fn="${2:-sine}"
n="${3:-64}"
out="output/${domain}_${source_fn}_n${n}.csv"
mkdir -p output
./compare2d "$domain" "$source_fn" "$n" "$out"
python3 plot_compare.py "$out"

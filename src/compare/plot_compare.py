#!/usr/bin/env python3
"""Plot a graphical comparison of the FDM and FEM solutions produced by
compare2d (see main.cpp).

Usage:
    python3 plot_compare.py result.csv [output.png]

Reads result.csv (node fields) and result.csv.residuals.csv (CG
convergence histories) and renders two figures:

  <output>.png      2x2: FDM field, FEM field, |difference|, CG convergence
  <output>_3d.png    1x2: FDM surface, FEM surface (both as 3D plots,
                     sharing one z-axis so their heights are comparable)
"""
import sys
import csv

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401  (registers the '3d' projection)


def read_header(path):
    with open(path) as fh:
        first = fh.readline()
    meta = {}
    if first.startswith("#"):
        for tok in first.lstrip("#").split():
            if "=" in tok:
                k, v = tok.split("=", 1)
                meta[k] = v
    return meta


def read_grid(path):
    with open(path) as fh:
        lines = [ln for ln in fh if not ln.startswith("#")]
    reader = csv.DictReader(lines)
    rows = list(reader)
    n = max(int(r["i"]) for r in rows)
    side = n + 1
    X = np.zeros((side, side))
    Y = np.zeros((side, side))
    active = np.zeros((side, side), dtype=bool)
    Ufdm = np.full((side, side), np.nan)
    Ufem = np.full((side, side), np.nan)
    for r in rows:
        i, j = int(r["i"]), int(r["j"])
        X[i, j] = float(r["x"])
        Y[i, j] = float(r["y"])
        act = r["active"] == "1"
        active[i, j] = act
        if act:
            Ufdm[i, j] = float(r["u_fdm"])
            Ufem[i, j] = float(r["u_fem"])
    return X, Y, active, Ufdm, Ufem


def plot_2d(meta, X, Y, Ufdm, Ufem, csv_path, out_path):
    diff = np.abs(Ufdm - Ufem)

    fig, axes = plt.subplots(2, 2, figsize=(10, 9))
    title = "Poisson -Laplacian(u) = f    domain={domain}  f={source}  N={n}".format(
        domain=meta.get("domain", "?"), source=meta.get("source", "?"), n=meta.get("n", "?")
    )
    fig.suptitle(title, fontsize=12)

    vmin = np.nanmin([np.nanmin(Ufdm), np.nanmin(Ufem)])
    vmax = np.nanmax([np.nanmax(Ufdm), np.nanmax(Ufem)])

    ax = axes[0, 0]
    im = ax.pcolormesh(X, Y, Ufdm, shading="auto", cmap="viridis", vmin=vmin, vmax=vmax)
    ax.set_title("Finite differences (5-point)")
    ax.set_aspect("equal")
    fig.colorbar(im, ax=ax)

    ax = axes[0, 1]
    im = ax.pcolormesh(X, Y, Ufem, shading="auto", cmap="viridis", vmin=vmin, vmax=vmax)
    ax.set_title("Finite elements (P1)")
    ax.set_aspect("equal")
    fig.colorbar(im, ax=ax)

    ax = axes[1, 0]
    im = ax.pcolormesh(X, Y, diff, shading="auto", cmap="magma")
    ax.set_title("|u_fdm - u_fem|  (max = {:.2e})".format(np.nanmax(diff)))
    ax.set_aspect("equal")
    fig.colorbar(im, ax=ax)

    ax = axes[1, 1]
    try:
        res_path = csv_path + ".residuals.csv"
        it, rfdm, rfem = [], [], []
        with open(res_path) as fh:
            reader = csv.DictReader(fh)
            for row in reader:
                it.append(int(row["iteration"]))
                rfdm.append(float(row["fdm_residual"]))
                rfem.append(float(row["fem_residual"]))
        ax.semilogy(it, rfdm, label="FDM")
        ax.semilogy(it, rfem, label="FEM")
        ax.set_xlabel("CG iteration")
        ax.set_ylabel("relative residual")
        ax.set_title("Conjugate-gradient convergence")
        ax.legend()
        ax.grid(True, which="both", alpha=0.3)
    except FileNotFoundError:
        ax.axis("off")

    fig.tight_layout(rect=[0, 0, 1, 0.95])
    fig.savefig(out_path, dpi=150)
    print("Wrote", out_path)


def plot_3d(meta, X, Y, Ufdm, Ufem, out_path):
    # NaN outside the domain gives a true gap in the surface instead of a
    # flat skirt at z=0, so the shape of the domain is still visible.
    vmin = np.nanmin([np.nanmin(Ufdm), np.nanmin(Ufem)])
    vmax = np.nanmax([np.nanmax(Ufdm), np.nanmax(Ufem)])

    fig = plt.figure(figsize=(11, 5.5))
    title = "Poisson -Laplacian(u) = f    domain={domain}  f={source}  N={n}".format(
        domain=meta.get("domain", "?"), source=meta.get("source", "?"), n=meta.get("n", "?")
    )
    fig.suptitle(title, fontsize=12)

    for idx, (label, Z) in enumerate([("Finite differences (5-point)", Ufdm), ("Finite elements (P1)", Ufem)]):
        ax = fig.add_subplot(1, 2, idx + 1, projection="3d")
        surf = ax.plot_surface(X, Y, Z, cmap="viridis", vmin=vmin, vmax=vmax,
                                linewidth=0, antialiased=True, rstride=1, cstride=1)
        ax.set_title(label)
        ax.set_xlabel("x")
        ax.set_ylabel("y")
        ax.set_zlabel("u")
        ax.set_zlim(vmin, vmax)
        ax.view_init(elev=28, azim=-60)
        fig.colorbar(surf, ax=ax, shrink=0.6, pad=0.08)

    fig.tight_layout(rect=[0, 0, 1, 0.95])
    fig.savefig(out_path, dpi=150)
    print("Wrote", out_path)


def main():
    if len(sys.argv) < 2:
        print(__doc__)
        return 1
    csv_path = sys.argv[1]
    out_path = sys.argv[2] if len(sys.argv) > 2 else csv_path.rsplit(".", 1)[0] + ".png"
    out_path_3d = out_path.rsplit(".", 1)[0] + "_3d.png"

    meta = read_header(csv_path)
    X, Y, active, Ufdm, Ufem = read_grid(csv_path)

    plot_2d(meta, X, Y, Ufdm, Ufem, csv_path, out_path)
    plot_3d(meta, X, Y, Ufdm, Ufem, out_path_3d)
    return 0


if __name__ == "__main__":
    sys.exit(main())

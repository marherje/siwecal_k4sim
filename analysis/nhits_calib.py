#!/usr/bin/env python
"""
Energy as a function of nhits, data and simulation: E = alpha n + beta n^2, from the per-point means written by
analysis/energy_trend.py (energy_trend.npz). Shows how much of a quadratic nhits->energy correction is plain pad
saturation (present in both) and how much differs between data and simulation.

    python analysis/nhits_calib.py <energy_trend dir>

Per set: power law nhits ~ E^p; least-squares E = alpha n + beta n^2 (and E = alpha n), parameter errors from
the fit residuals (unit weights, the beam energies taken as exact). Writes nhits_calib.txt and nhits_calib.png.
"""
import os
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

DATA_C, SIM_C = "#0b0b0b", "#2a78d6"


def fit(n, E, cols):
    A = np.vstack([c(n) for c in cols]).T
    p, *_ = np.linalg.lstsq(A, E, rcond=None)
    r = A @ p - E
    dof = max(len(E) - len(p), 1)
    s2 = (r ** 2).sum() / dof
    pe = np.sqrt(np.diag(np.linalg.inv(A.T @ A)) * s2)
    return p, pe, np.sqrt(np.mean(r ** 2))


def main():
    d = sys.argv[1]
    z = np.load(f"{d}/energy_trend.npz")
    sets = sorted({k.split("|")[0] for k in z.files})
    lines = []
    plt.rcParams.update({"axes.spines.top": False, "axes.spines.right": False, "axes.grid": True,
                         "grid.color": "#e6e6e3", "axes.edgecolor": "#8c8c8c", "font.size": 10})
    fig, ax = plt.subplots(2, len(sets), figsize=(5.6 * len(sets), 6.6), squeeze=False,
                           gridspec_kw={"height_ratios": [2.2, 1]}, sharex="col")
    for j, s in enumerate(sets):
        E = z[f"{s}|E"]
        lines.append(f"== {s}   E: {' '.join(f'{e:g}' for e in E)}")
        for lab, key, col, mk in (("data", "nh_d", DATA_C, "o"), ("simulation", "nh_s", SIM_C, "s")):
            n = z[f"{s}|{key}"]
            pw = np.polyfit(np.log(E), np.log(n), 1)[0]
            (a2, b2), (ea2, eb2), rms2 = fit(n, E, [lambda x: x, lambda x: x ** 2])
            (a1,), (ea1,), rms1 = fit(n, E, [lambda x: x])
            share = b2 * n.max() ** 2 / (a2 * n.max() + b2 * n.max() ** 2)
            lines.append(f"  {lab:<10} nhits ~ E^{pw:.3f} | E = ({a2:.4f} +- {ea2:.4f}) n + ({b2 * 1e4:.3f} +- {eb2 * 1e4:.3f})e-4 n^2"
                         f"  rms {rms2:.2f} GeV | E = {a1:.4f} n  rms {rms1:.2f} GeV | n^2 share at max n {share * 100:.0f} %")
            nf = np.linspace(0, n.max() * 1.05, 200)
            ax[0, j].plot(n, E, mk, ms=7, color=col, label=f"{lab}  (nhits ∝ E^{pw:.2f})", zorder=3)
            ax[0, j].plot(nf, a2 * nf + b2 * nf ** 2, "-", lw=2, color=col,
                          label=f"  E = {a2:.4f}·n + {b2 * 1e4:.2f}e-4·n²")
            ax[0, j].plot(nf, a1 * nf, ":", lw=1.5, color=col)
            ax[1, j].plot(n, (a2 * n + b2 * n ** 2) - E, mk + "-", ms=6, lw=1.5, color=col)
        ax[0, j].set_title(s + " (barycentre > 30 mm from edges)", fontsize=10, loc="left")
        ax[0, j].set_ylabel("beam energy [GeV]"); ax[0, j].legend(frameon=False, fontsize=8, loc="upper left")
        ax[1, j].axhline(0, color="#52514e", lw=1); ax[1, j].set_ylabel("fit − E [GeV]")
        ax[1, j].set_xlabel("hits per event (mean)")
        bd = [float(l.split("(")[2].split("+-")[0]) for l in lines[-2:]]
        lines.append(f"  beta(data) / beta(sim) = {bd[0] / bd[1]:.2f}")
    fig.text(0.01, 0.005, "solid: E = α·n + β·n²   dotted: E = α·n", fontsize=8, color="#52514e")
    fig.tight_layout(); fig.savefig(f"{d}/nhits_calib.png", dpi=130); plt.close(fig)
    open(f"{d}/nhits_calib.txt", "w").write("\n".join(lines) + "\n")
    print("\n".join(lines))


if __name__ == "__main__":
    main()

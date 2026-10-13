#!/usr/bin/env python
"""
Energy linearity, data against the digitised simulation (final_v6): Gaussian-core mean of the shower energy against
the beam energy, the deviation from a straight line fitted at low energy, and the sim / data ratio per energy.

    FINAL_SIM=.../final_v6 python analysis/linearity.py <outdir> [--cache-dir analysis/out_panels_v6]

Samples: th230 E4 scan 7.5-74 GeV (data runs 41-46), and th210 S1 at 20 and 74 GeV (eudaq 291, 286).
Inputs and selection are final_panels.py's (valtrees, 0.5 MIP hit cut, shower-like events); the per-group cache of
final_panels.py is reused. The straight line E_reco = a E + b is fitted per sample kind on 7.5-20 GeV, so the
residual at higher energy is the non-linearity of each kind on its own scale.

Writes <outdir>/linearity.png and linearity.txt.
"""
import argparse
import os
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(_HERE))
os.environ.setdefault("FINAL_GROUPS", os.path.join(_HERE, "final_v6_groups.txt"))
from analysis.final_panels import COLOR, LABEL, group_samples, load_groups  # noqa: E402
from analysis.compare_event_level import gauss_core  # noqa: E402

INK, MUTED, GRID = "#1b2227", "#56626b", "#d7dde2"
plt.rcParams.update({"font.size": 11, "axes.edgecolor": MUTED, "axes.labelcolor": INK, "xtick.color": MUTED,
                     "ytick.color": MUTED, "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.7,
                     "axes.spines.top": False, "axes.spines.right": False, "legend.frameon": False})
SCANS = [("th230 E4", "th230_E4_e", ["7.5", "10", "20", "34", "52", "74"]), ("th210 S1", "th210_S1_e", ["20", "74"])]
VARS = [("sum_energy", "Energy [MIP]"), ("sum_w", "W-weighted energy [MIP]"), ("nhit_chan", "Number of hits")]
KINDS = ["data", "digi"]
FIT_MAX = 20.0


def measure(cache_dir):
    groups = {g["name"]: g for g in load_groups(os.environ["FINAL_GROUPS"])}
    out = {}
    for scan, pre, energies in SCANS:
        for E in energies:
            s = group_samples(groups[pre + E], cache_dir)
            for k in ("data", "nodigi", "digi"):
                for v, _ in VARS:
                    x = s[k][v][np.isfinite(s[k][v])]
                    mu, sg = gauss_core(x)
                    out[(scan, float(E), k, v)] = (mu, sg, sg / np.sqrt(max(len(x), 1)), len(x))
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("outdir")
    ap.add_argument("--cache-dir", default=os.path.join(_HERE, "out_panels_v6"))
    a = ap.parse_args()
    os.makedirs(a.outdir, exist_ok=True)
    m = measure(a.cache_dir)
    fig, ax = plt.subplots(3, len(VARS), figsize=(5.2 * len(VARS), 12), sharex=True,
                           gridspec_kw={"height_ratios": [2, 1, 1]})
    lines = [f"Gaussian-core mean per sample; line a E + b fitted on 7.5-{FIT_MAX:.0f} GeV per kind (th230 E4)"]
    for j, (v, yl) in enumerate(VARS):
        lines.append(f"\n== {v}")
        for scan, _, energies in SCANS:
            E = np.array([float(e) for e in energies])
            mk = "o" if scan == "th230 E4" else "D"
            mu = {k: np.array([m[(scan, e, k, v)][0] for e in E]) for k in KINDS}
            er = {k: np.array([m[(scan, e, k, v)][2] for e in E]) for k in KINDS}
            for k in KINDS:
                fill = COLOR[k] if scan == "th230 E4" else "none"
                ax[0, j].errorbar(E, mu[k], yerr=er[k], fmt=mk, color=COLOR[k], mfc=fill, ms=7,
                                  label=f"{LABEL[k]}, {scan}")
                if scan == "th230 E4":
                    lo = E <= FIT_MAX
                    p = np.polyfit(E[lo], mu[k][lo], 1)
                    xx = np.linspace(0, 80, 50)
                    ax[0, j].plot(xx, np.polyval(p, xx), "--", color=COLOR[k], lw=1.2)
                    res = mu[k] / np.polyval(p, E) - 1
                    ax[1, j].errorbar(E, 100 * res, yerr=100 * er[k] / np.polyval(p, E), fmt=mk + "-",
                                      color=COLOR[k], ms=6, label=LABEL[k])
                    lines.append(f"  {LABEL[k]:<18} a = {p[0]:8.2f}/GeV  b = {p[1]:8.1f}   residual [%]: "
                                 + " ".join(f"{e:g}:{100 * r:+.1f}" for e, r in zip(E, res)))
            r = mu["digi"] / mu["data"]
            rer = r * np.hypot(er["digi"] / mu["digi"], er["data"] / mu["data"])
            ax[2, j].errorbar(E, r, yerr=rer, fmt=mk + ("-" if scan == "th230 E4" else ""), color=COLOR["digi"],
                              mfc=COLOR["digi"] if scan == "th230 E4" else "none", ms=6, label=scan)
            lines.append(f"  {scan}: sim/data " + " ".join(f"{e:g}:{x:.3f}" for e, x in zip(E, r))
                         + "   data mu: " + " ".join(f"{x:.0f}" for x in mu["data"])
                         + "   sim mu: " + " ".join(f"{x:.0f}" for x in mu["digi"]))
        ax[0, j].set_ylabel(yl); ax[0, j].set_title(yl, loc="left")
        ax[1, j].set_ylabel(f"deviation from line (fit $\\leq${FIT_MAX:.0f} GeV) [%]"); ax[1, j].axhline(0, color=INK, lw=0.8)
        ax[2, j].set_ylabel("sim / data"); ax[2, j].axhline(1, color=INK, lw=0.8)
        ax[2, j].set_xlabel("beam energy [GeV]")
        ax[1, j].legend(fontsize=8.5); ax[2, j].legend(fontsize=8.5)
    ax[0, 0].legend(fontsize=8.5, loc="upper left")
    for x in ax.ravel():
        x.set_xlim(0, 80)
    fig.tight_layout(); fig.savefig(f"{a.outdir}/linearity.png", dpi=120); plt.close(fig)
    txt = "\n".join(lines)
    print(txt)
    open(f"{a.outdir}/linearity.txt", "w").write(txt + "\n")


if __name__ == "__main__":
    main()

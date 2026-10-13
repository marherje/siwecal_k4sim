#!/usr/bin/env python
"""
Figures of analysis/cm_trigger_test.py: hit-bit turn-on by chip-SCA charge, in the raw and the baseline-corrected
residual, and the fitted 50 % point against the charge.

    python analysis/cm_trigger_plot.py <dir> <name> <label>=<tag> [<label>=<tag> ...]

<tag> is the suffix of cm_trigger_<tag>.npz (e291, e286_6_11, ...). Writes <dir>/<name>.png and <name>.txt.
"""
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

from cm_trigger_test import fit

INK, MUTED, GRID = "#1b2227", "#56626b", "#d7dde2"
RAMP = ["#c6dbef", "#9ecae1", "#6baed6", "#4292c6", "#2171b5", "#08519c", "#08306b"]
plt.rcParams.update({"font.size": 11, "axes.edgecolor": MUTED, "axes.labelcolor": INK, "xtick.color": MUTED,
                     "ytick.color": MUTED, "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.7,
                     "axes.spines.top": False, "axes.spines.right": False, "legend.frameon": False})


def main():
    d, name, items = sys.argv[1], sys.argv[2], [x.split("=") for x in sys.argv[3:]]
    fig, ax = plt.subplots(len(items), 3, figsize=(15, 4.2 * len(items)), squeeze=False)
    lines = []
    for row, (lab, tag) in enumerate(items):
        z = np.load(f"{d}/cm_trigger_{tag}.npz")
        qe, rb, base = z["q_edges"], z["r_bins"], z["base"]
        ctr = 0.5 * (rb[1:] + rb[:-1])
        qc, mr, mc, er = [], [], [], []
        lines.append(f"== {lab} ({tag}, {int(z['n_acq'])} acquisitions, {int(z['masked'].sum())} channels excluded)")
        lines.append(f"{'Q bin':>13} {'<B>':>6} {'raw mu':>7} {'cor mu':>7} {'plateau':>8}")
        for k in range(len(qe) - 1):
            if z["raw_all"][k][(ctr > 5) & (ctr < 35)].sum() < 1500:   # too few samples in the turn-on region
                continue
            fr = fit(z["raw_bit"][k], z["raw_all"][k])
            fc = fit(z["cor_bit"][k], z["cor_all"][k])
            if not np.isfinite(fr[2]):
                continue
            lo, hi = qe[k], qe[k + 1]
            qlab = f"{lo:.0f}-{hi:.0f}" if hi < 1e8 else f">{lo:.0f}"
            lines.append(f"{qlab:>13} {base[k]:+6.2f} {fr[2]:7.2f} {fc[2]:7.2f} {fr[1]:8.2f}")
            col = RAMP[min(k, len(RAMP) - 1)]
            for j, (key, f) in enumerate((("raw", fr), ("cor", fc))):
                a, b = z[f"{key}_all"][k], z[f"{key}_bit"][k]
                ok = a > 30
                ax[row, j].plot(ctr[ok], (b[ok] / a[ok] - f[0]) / (f[1] - f[0]), "-", color=col, lw=1.6,
                                label=f"Q {qlab}, B {base[k]:+.1f}")
            qc.append(np.sqrt(max(lo, 50) * min(hi, 2 * max(lo, 50)))); mr.append(fr[2]); mc.append(fc[2])
            er.append(base[k])
        qc, mr, mc, er = map(np.array, (qc, mr, mc, er))
        ax[row, 2].plot(qc, mr - mr[0], "o-", color="#2a78d6", label="50 % point, raw ADC - pedestal")
        ax[row, 2].plot(qc, mc - mc[0], "s-", color="#c0392b", label="50 % point, corrected for the baseline")
        ax[row, 2].plot(qc, er - er[0], "^:", color=MUTED, label="baseline B")
        ax[row, 2].axhline(0, color=INK, lw=0.8)
        ax[row, 2].set_xscale("log")
        ax[row, 2].set_xlabel("charge of the other fired channels in the chip-SCA [HG ADC]")
        ax[row, 2].set_ylabel("shift from the quietest bin [ADC]")
        ax[row, 2].legend(fontsize=8.5)
        for j, xl in enumerate(("ADC - pedestal (raw)", "ADC - pedestal - baseline (corrected)")):
            ax[row, j].set_xlim(0, 40); ax[row, j].set_ylim(-0.05, 1.1)
            ax[row, j].set_xlabel(xl); ax[row, j].set_ylabel("hit-bit efficiency / plateau")
            ax[row, j].axhline(0.5, color=INK, lw=0.6, ls=":")
        ax[row, 0].legend(fontsize=8, loc="lower right")
        ax[row, 0].set_title(lab, loc="left", fontsize=11.5)
    fig.tight_layout(); fig.savefig(f"{d}/{name}.png", dpi=120); plt.close(fig)
    txt = "\n".join(lines)
    print(txt)
    open(f"{d}/{name}.txt", "w").write(txt + "\n")


if __name__ == "__main__":
    main()

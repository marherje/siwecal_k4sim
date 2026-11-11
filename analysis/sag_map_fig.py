#!/usr/bin/env python
"""
Figure of the measured baseline-sag map against the data, per threshold set (upstream-line note, section 2.4).

Top row: mean baseline B of signal-free channels against Q_chip, in four slices of Q_slab; data as points (one marker
per slice), the measured map (occupancy_sag_<set>.txt, evaluated on the same records) as lines (one line style per
slice). Bottom row, one panel: data - map (dark red, markers) and data - k*Q_chip, the previous digitiser
(k = -1.5e-3; light red, thin lines), one curve per slice with the slice's line style.

    python analysis/sag_map_fig.py <dir with cm_sag_*.npz> <dir with occupancy_sag_<set>.txt> <out.png>
"""
import argparse
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from analysis import plot_style as ps  # noqa: E402
from analysis.sag_model_fit import load, judge  # noqa: E402
from analysis.sag_table import SETS, QC_EDGES, interp  # noqa: E402

K_PREV = -1.5e-3
SLICES = [(0, 2000, "Q_slab < 2k"), (2000, 8000, "Q_slab 2–8k"), (8000, 16000, "Q_slab 8–16k"),
          (16000, np.inf, "Q_slab > 16k")]
MARK = ["o", "s", "^", "D"]
DATA_COL = ["#29a3dc", "#1f7fb0", "#7cc6ea", "#0d4a6b"]
MIN_N = 50
ADC_PER_MIP = 18.75
ANALYTIC = "#7a7a7a"   # grey: a formula, neither data nor the digitiser's model
TITLE = {"th210": "th210 (20, 52, 74 GeV)", "th220": "th220 (20, 74 GeV)", "th230": "th230 (20, 52, 74 GeV)"}


def read_table(path):
    xc = yc = None
    rows = []
    for line in open(path):
        if line.startswith("#"):
            continue
        f = line.split()
        if f[0] == "qchip":
            xc = np.array(f[1:], float)
        elif f[0] == "qslab":
            yc = np.array(f[1:], float)
        else:
            rows.append([float(v) for v in f])
    return xc, yc, np.array(rows)


Q0_GRID = [10, 30, 100, 300, 1000, 3000]


def fit_analytic(qc, qs, B):
    """Best analytic model of the study: B = (a + e Q_slab) ln(1 + Q_chip/q0) + b Q_slab + k Q_chip, least squares on
    the records (linear in a, e, b, k; q0 from a grid). Returns (prediction function, parameters)."""
    best = None
    for q0 in Q0_GRID:
        L = np.log1p(qc / q0)
        X = np.c_[L, qs * L, qs, qc]
        p = np.linalg.lstsq(X, B, rcond=None)[0]
        rss = ((B - X @ p) ** 2).sum()
        if best is None or rss < best[0]:
            best = (rss, q0, p)
    _, q0, (a_, e_, b_, k_) = best

    def f(qc_, qs_):
        L = np.log1p(qc_ / q0)
        return (a_ + e_ * qs_) * L + b_ * qs_ + k_ * qc_
    return f, dict(a=a_, e=e_, b=b_, k=k_, q0=q0)


def profile(qc, y):
    i = np.digitize(qc, QC_EDGES) - 1
    out = []
    for b in range(len(QC_EDGES) - 1):
        m = i == b
        if m.sum() >= MIN_N:
            out.append((qc[m].mean(), y[m].mean(), y[m].std() / np.sqrt(m.sum())))
    return np.array(out).reshape(-1, 3)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("records")
    ap.add_argument("tables")
    ap.add_argument("out")
    a = ap.parse_args()
    c = load(a.records)
    plt.rcParams.update(ps.RC)
    fig, ax = plt.subplots(2, 3, figsize=(17, 9.5), sharex=True, sharey="row",
                           gridspec_kw=dict(height_ratios=[2.4, 1.3]))
    for col, (th, runs) in enumerate(SETS.items()):
        m = np.isin(c["run"], runs)
        qc, qs, B = c["q_chip"][m], c["q_slab"][m], c["base_far"][m]
        xc, yc, tab = read_table(f"{a.tables}/occupancy_sag_{th}.txt")
        pred_map = interp(xc, yc, tab, qc, qs)
        pred_kq = K_PREV * qc
        fa, par = fit_analytic(qc, qs, B)
        pred_an = fa(qc, qs)
        cc = {k: v[m] for k, v in c.items()}
        X2 = np.c_[qc, qs]
        pred_lin = X2 @ np.linalg.lstsq(X2, B, rcond=None)[0]
        print(f"{th}: best formula a={par['a']:.3f} e={par['e']:.3e} b={par['b']:.3e} k={par['k']:.3e} q0={par['q0']:g}")
        print(f"{th}: chi2/ndf (sag_model_fit.judge bins, per run): kQ {judge(cc, pred_kq)[1]:.0f}, "
              f"kQ+bQs {judge(cc, pred_lin)[1]:.0f}, formula {judge(cc, pred_an)[1]:.0f}, map {judge(cc, pred_map)[1]:.0f}")
        for s, (lo, hi, lab) in enumerate(SLICES):
            sel = (qs >= lo) & (qs < hi)
            if sel.sum() < MIN_N:
                continue
            d = profile(qc[sel], B[sel])
            pm = profile(qc[sel], pred_map[sel])
            pa = profile(qc[sel], pred_an[sel])
            ra = profile(qc[sel], B[sel] - pred_an[sel])
            rm = profile(qc[sel], B[sel] - pred_map[sel])
            rk = profile(qc[sel], B[sel] - pred_kq[sel])
            ls = ps.LINESTYLES[s]
            ax[0, col].errorbar(d[:, 0], d[:, 1], yerr=d[:, 2], fmt=MARK[s], color=DATA_COL[s], ms=5, mew=0,
                                elinewidth=1.2, label=f"data, {lab}", zorder=5)
            ax[0, col].plot(pm[:, 0], pm[:, 1], color=ps.DIGI, ls=ls, lw=1.6, label=f"measured map, {lab}")
            ax[0, col].plot(pa[:, 0], pa[:, 1], color=ANALYTIC, ls=ls, lw=1.3, label=f"best formula, {lab}")

            ax[1, col].errorbar(rm[:, 0], rm[:, 1], yerr=rm[:, 2], fmt=MARK[s], color=ps.DIGI, ms=4, mew=0,
                                elinewidth=1, ls=ls, lw=1.2)
            ax[1, col].plot(rk[:, 0], rk[:, 1], color=ps.REDS[0], ls=ls, lw=1.1, zorder=1)
            ax[1, col].plot(ra[:, 0], ra[:, 1], color=ANALYTIC, ls=ls, lw=1.3, zorder=2)
        qq = np.logspace(0, np.log10(1.5e4), 100)
        ax[0, col].plot(qq, K_PREV * qq, color=ps.REDS[0], lw=1.3, zorder=1,
                        label="k·Q_chip, k = −1.5·10⁻³\n(same for every Q_slab)")
        ax[0, col].set_title(TITLE[th])
        ax[1, col].axhline(0, color="0.3", lw=0.8)
        ax[1, col].set_ylim(-7, 7)
        # reference lines in MIP (AdcPerMip = 18.75): top -0.25, -0.5 and -0.75 MIP, bottom +0.25 and -0.25 MIP
        for v, lab in ((-0.25, "−0.25 MIP"), (-0.5, "−0.5 MIP"), (-0.75, "−0.75 MIP")):
            ax[0, col].axhline(v * ADC_PER_MIP, color="0.45", lw=0.9, ls=(0, (6, 3)))
            ax[0, col].text(1.15, v * ADC_PER_MIP + 0.25, lab, fontsize=8, color="0.35")
        for v in (-0.25, 0.25):
            ax[1, col].axhline(v * ADC_PER_MIP, color="0.45", lw=0.9, ls=(0, (6, 3)))
        ax[1, col].text(1.8e4, 0.25 * ADC_PER_MIP + 0.3, "+0.25 MIP", fontsize=8, color="0.35", ha="right")
        ax[1, col].text(1.8e4, -0.25 * ADC_PER_MIP - 1.0, "−0.25 MIP", fontsize=8, color="0.35", ha="right")
        ax[1, col].set_xlabel("Q_chip, charge of the chip-SCA [ADC]")
        for r in range(2):
            ax[r, col].set_xscale("log")
            ax[r, col].set_xlim(1, 2e4)
            ax[r, col].grid(alpha=0.25)
    ax[0, 0].set_ylabel("Δ Pedestal (ADC)")
    ax[0, 0].set_ylim(-15, 1)
    ax[1, 0].set_ylabel("data − model [ADC]")
    from matplotlib.lines import Line2D
    ax[1, 0].legend(handles=[Line2D([], [], color=ps.DIGI, marker="o", ms=4, lw=1.2, label="data − measured map"),
                             Line2D([], [], color=ANALYTIC, lw=1.3, label="data − best formula"),
                             Line2D([], [], color=ps.REDS[0], lw=1.1,
                                    label="data − k·Q_chip (k = −1.5·10⁻³)")],
                    loc="upper left", fontsize=8, frameon=True, framealpha=0.9)
    h, l = ax[0, 0].get_legend_handles_labels()
    dat = [i for i in range(len(l)) if l[i].startswith("data")]
    mod = [i for i in range(len(l)) if l[i].startswith("meas")]
    frm = [i for i in range(len(l)) if l[i].startswith("best")]
    kq = [i for i in range(len(l)) if l[i].startswith("k·Q")]
    order = [x for tri in zip(dat, mod, frm) for x in tri]   # column-major: rows data, map, formula
    fig.legend([h[i] for i in order], [l[i] for i in order], loc="upper center", bbox_to_anchor=(0.45, 0.975), ncol=4,
               frameon=False, fontsize=9)
    # k*Q, one line for every slice: its own legend, set apart on the right
    fig.legend([h[kq[0]]], [l[kq[0]]], loc="upper left", bbox_to_anchor=(0.738, 0.965), frameon=True, fontsize=9)
    fig.suptitle("Chip baseline sag against the charge of the chip and of the rest of its slab, per threshold set",
                 y=0.995, fontsize=12)
    fig.tight_layout(rect=[0, 0, 1, 0.87])
    fig.savefig(a.out, dpi=110)
    print(a.out)


if __name__ == "__main__":
    main()

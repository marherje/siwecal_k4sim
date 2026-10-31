#!/usr/bin/env python
"""
sim_tuning 2.3: does the high-gain roll-over, or the LG/HG anchor, change with the charge in the chip?

For every hit with HG > 150 ADC: the ratio LG/HG against HG, in bins of the chip charge Q (energy of the other hits
of the same chip in the event, MIP), data and digitised simulation of the same beam points. A gain that falls with
Q shows as a lower plateau (anchor); a roll-over that sets in earlier shows as the ratio rising at lower HG.

Per Q bin: plateau = median LG/HG at HG 150-800 ADC; roll-over point A = HG where the median ratio exceeds the
plateau by 5 % (the high gain falling behind). Writes <out>/gain_<group>.png, gain.txt and gain_rq.npz (plateau
against Q for apply_lowgain).

    python analysis/tuning_gain.py <outdir> [groups]
"""
import os
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from analysis.sim_tuning import load, load_groups, data_paths, sim_paths, chip_charge, OUT  # noqa: E402

Q_EDGES = np.array([0, 5, 15, 30, 60, 120, 250, 500, 5000.0])
HG_EDGES = np.array([150, 200, 300, 400, 550, 700, 850, 1000, 1150, 1300, 1450, 1600, 1750, 1900, 2050, 2200.0])
COL = plt.cm.viridis(np.linspace(0.05, 0.9, len(Q_EDGES) - 1))


def ratio_table(s):
    q = chip_charge(s)
    hg, lg = s["hg"].astype(float), s["lg"].astype(float)
    k = hg > 150
    qb = np.digitize(q[k], Q_EDGES) - 1
    hb = np.digitize(hg[k], HG_EDGES) - 1
    r = lg[k] / hg[k]
    med = np.full((len(Q_EDGES) - 1, len(HG_EDGES) - 1), np.nan)
    n = np.zeros_like(med)
    for i in range(len(Q_EDGES) - 1):
        for j in range(len(HG_EDGES) - 1):
            m = (qb == i) & (hb == j)
            n[i, j] = m.sum()
            if m.sum() >= 30:
                med[i, j] = np.median(r[m])
    plat = np.array([np.median(r[(qb == i) & (hg[k] < 800)]) if ((qb == i) & (hg[k] < 800)).sum() > 30 else np.nan
                     for i in range(len(Q_EDGES) - 1)])
    hc = 0.5 * (HG_EDGES[1:] + HG_EDGES[:-1])
    A = []
    for i in range(len(Q_EDGES) - 1):
        above = np.nonzero(med[i] > 1.05 * plat[i])[0]
        A.append(hc[above[0]] if len(above) else np.nan)
    return med, n, plat, np.array(A)


def main():
    out = sys.argv[1]
    want = sys.argv[2:] or ["th210_e20_26_-52", "th210_e74_22_-49", "th230_E4_e20", "th230_E4_e74"]
    os.makedirs(out, exist_ok=True)
    groups = {g["name"]: g for g in load_groups()}
    qc = np.sqrt(np.maximum(Q_EDGES[:-1], 1) * Q_EDGES[1:])
    lines = ["group sample  plateau(Q bins)  |  roll-over A [HG ADC] (Q bins);  Q edges " + " ".join(map(str, Q_EDGES))]
    store = {}
    for n in want:
        g = groups[n]
        d = load(data_paths(g), cache=f"{OUT}/cache/data_{n}.npz")
        s = load(sim_paths(n), cache=f"{OUT}/cache/v7_{n}.npz")
        fig, ax = plt.subplots(1, 3, figsize=(16, 4.6))
        for lab, smp, a in (("data", d, ax[0]), ("sim", s, ax[1])):
            med, cnt, plat, A = ratio_table(smp)
            store[f"{n}_{lab}"] = plat
            hc = 0.5 * (HG_EDGES[1:] + HG_EDGES[:-1])
            for i in range(len(Q_EDGES) - 1):
                a.plot(hc, med[i], "-o", ms=3, color=COL[i], label=f"Q {Q_EDGES[i]:g}-{Q_EDGES[i + 1]:g} MIP")
            a.set_title(f"{n}: {lab}"); a.set_xlabel("HG [ADC]"); a.set_ylabel("median LG / HG")
            ax[2].plot(qc, plat, "-o" if lab == "data" else "--s", label=f"{lab} plateau (HG 150-800)")
            lines.append(f"{n} {lab}  " + " ".join(f"{x:.4f}" for x in plat) + "  |  " +
                         " ".join(f"{x:.0f}" for x in A))
        ax[0].legend(fontsize=7)
        ax[2].set_xscale("log"); ax[2].set_xlabel("chip charge Q of the other hits [MIP]")
        ax[2].set_ylabel("LG/HG plateau"); ax[2].legend(fontsize=8)
        fig.tight_layout(); fig.savefig(f"{out}/gain_{n}.png", dpi=110); plt.close(fig)
        print("\n".join(lines[-2:]), flush=True)
    open(f"{out}/gain.txt", "w").write("\n".join(lines) + "\n")
    np.savez(f"{out}/gain_rq.npz", q_c=qc, **store)


if __name__ == "__main__":
    main()

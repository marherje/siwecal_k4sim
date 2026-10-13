#!/usr/bin/env python
"""
Hit loss of every chip against its occupancy: one panel per slab (ASU), one curve per chip. Shower-like events of one
beam group from the valtrees (0.5 MIP hit cut in data and simulation).

The occupancy measured in data is itself reduced by the loss, so the x axis is the occupancy the SIMULATION gives at
the same quantile of the per-chip hit-count distribution (quantile matching): with the events of each sample sorted by
the number of hits the chip has, the event at quantile q has n_sim(q) hits in simulation and n_data(q) in data, and

    loss(q) = 1 - n_data(q) / n_sim(q)      plotted against n_sim(q).

This assumes the same population of showers in both samples (same beam spot and size), so any sim/data difference in
the lateral shower shape or beam position is read as loss too ("apparent" loss). Counts are smoothed with a uniform
[0, 1) jitter before taking quantiles; the band is the bootstrap spread (50 resamplings).

    python analysis/chip_loss.py --out <dir> --label "th210 S1 20 GeV" --data <data valtree> --sim <sim valtree>

Writes chip_loss_<tag>.png (3 x 5 panels) and chip_loss_<tag>.txt.
"""
import argparse

import numpy as np
import uproot
import awkward as ak
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402

from trigger_maps import showerlike, C_GRID, N_SLAB, N_CHIP

# colour = chip column (x), line style = chip row (y), both ranked within the slab
C_COL = ["#2a78d6", "#1baf7a", "#eda100", "#c0392b"]
LS_ROW = ["-", "--", "-.", ":"]
MK_ROW = ["o", "s", "^", "D"]
QS = np.linspace(0.40, 0.995, 36)
MIN_NSIM = 2.0          # quantiles where the simulation has fewer hits in the chip are not drawn


def load(path, nmax=0):
    a = uproot.open(path)["ecal"].arrays(["nhit_chan", "hit_slab", "hit_chip", "hit_x", "hit_y"],
                                         entry_stop=nmax or None, library="ak")
    a = a[showerlike(ak.to_numpy(a["nhit_chan"]).astype(float))]
    key = a["hit_slab"] * N_CHIP + a["hit_chip"]
    n = np.zeros((len(a), N_SLAB * N_CHIP))
    ev = np.repeat(np.arange(len(a)), ak.to_numpy(ak.num(key)))
    k = ak.to_numpy(ak.flatten(key)).astype(int)
    np.add.at(n, (ev, k), 1)
    x = np.bincount(k, ak.to_numpy(ak.flatten(a["hit_x"])), minlength=N_SLAB * N_CHIP)
    y = np.bincount(k, ak.to_numpy(ak.flatten(a["hit_y"])), minlength=N_SLAB * N_CHIP)
    c = np.bincount(k, minlength=N_SLAB * N_CHIP)
    with np.errstate(invalid="ignore", divide="ignore"):
        pos = np.stack([x / c, y / c], axis=1)
    return n, pos


def loss_curve(nd, ns, rng):
    jd = nd + rng.random(len(nd)); js = ns + rng.random(len(ns))
    qd, qs_ = np.quantile(jd, QS), np.quantile(js, QS)
    with np.errstate(invalid="ignore", divide="ignore"):
        return qs_, 1 - qd / qs_


def rank4(v):
    """Rank four-ish distinct coordinates of the chips of a slab into 0..3 (by the 4 clusters of chip centres)."""
    out = np.full(len(v), -1)
    ok = np.isfinite(v)
    if ok.sum() == 0:
        return out
    order = np.argsort(v[ok])
    bins = np.array_split(order, 4)
    idx = np.flatnonzero(ok)
    for r, b in enumerate(bins):
        out[idx[b]] = r
    return out


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", required=True)
    ap.add_argument("--label", required=True)
    ap.add_argument("--data", required=True)
    ap.add_argument("--sim", required=True)
    ap.add_argument("--sim-label", default="per-slab trigger model")
    ap.add_argument("--max-events", type=int, default=0)
    a = ap.parse_args(argv)
    tag = a.label.replace(" ", "_")
    nd, pd = load(a.data, a.max_events)
    ns, ps = load(a.sim, a.max_events)
    pos = np.where(np.isfinite(pd), pd, ps)
    rng = np.random.default_rng(1)

    fig, axes = plt.subplots(3, 5, figsize=(25, 14.5), sharex=True, sharey=True, constrained_layout=True)
    lines = [f"# {a.label}: {len(nd)} data, {len(ns)} sim shower-like events; loss [%] at n_sim quantiles {list(np.round(QS, 3))}"]
    xmax = 1
    for slab, ax in enumerate(axes.ravel()):
        sl = slice(slab * N_CHIP, (slab + 1) * N_CHIP)
        col = rank4(pos[sl, 0]); row = rank4(-pos[sl, 1])     # column left->right, row top->bottom
        drawn = 0
        for chip in range(N_CHIP):
            k = slab * N_CHIP + chip
            if ns[:, k].mean() < 0.05 or nd[:, k].sum() < 50:
                continue
            x, l = loss_curve(nd[:, k], ns[:, k], rng)
            # bootstrap band
            boots = []
            for _ in range(50):
                i = rng.integers(0, len(nd), len(nd)); j = rng.integers(0, len(ns), len(ns))
                boots.append(loss_curve(nd[i, k], ns[j, k], rng)[1])
            e = np.nanstd(np.array(boots), axis=0)
            m = x >= MIN_NSIM
            if m.sum() < 2:
                continue
            c = C_COL[col[chip]] if col[chip] >= 0 else "0.4"
            r = row[chip] if row[chip] >= 0 else 0
            ax.plot(x[m], 100 * l[m], LS_ROW[r], color=c, lw=1.8)
            ax.fill_between(x[m], 100 * (l[m] - e[m]), 100 * (l[m] + e[m]), color=c, alpha=0.12, lw=0)
            ax.plot(x[m][-1], 100 * l[m][-1], MK_ROW[r], color=c, ms=6)
            ax.annotate(str(chip), (x[m][-1], 100 * l[m][-1]), xytext=(4, 0), textcoords="offset points",
                        fontsize=9, color="0.2", va="center")
            xmax = max(xmax, x[m].max())
            drawn += 1
            lines.append(f"slab {slab:2d} chip {chip:2d} (x {pos[k, 0]:6.1f}, y {pos[k, 1]:6.1f}): "
                         + " ".join(f"{xx:.1f}:{100 * ll:+.1f}" for xx, ll in zip(x[m], l[m])))
        ax.axhline(0, color="0.45", lw=1)
        ax.grid(color=C_GRID, lw=0.6)
        ax.set_title(f"ASU / slab {slab}  ({drawn} chips)", fontsize=13)
    for ax in axes[-1]:
        ax.set_xlabel("hits in the chip in simulation (same quantile)", fontsize=12)
    for ax in axes[:, 0]:
        ax.set_ylabel("hits lost in data [%]", fontsize=12)
    axes[0, 0].set_xlim(MIN_NSIM - 0.5, xmax * 1.12)
    axes[0, 0].set_ylim(-40, 60)
    handles = [Line2D([], [], color=C_COL[i], lw=2.4, label=f"chip column {i + 1} (left to right)") for i in range(4)] + \
              [Line2D([], [], color="0.25", ls=LS_ROW[i], marker=MK_ROW[i], lw=1.8, label=f"chip row {i + 1} (top to bottom)")
               for i in range(4)]
    fig.legend(handles=handles, loc="outside lower center", ncol=8, fontsize=11, frameon=False)
    fig.suptitle(f"{a.label}: hits lost in data vs chip occupancy, quantile-matched to the simulation ({a.sim_label}); "
                 "number at the end of each curve = chip", fontsize=15)
    fig.savefig(f"{a.out}/chip_loss_{tag}.png", dpi=105); plt.close(fig)
    with open(f"{a.out}/chip_loss_{tag}.txt", "w") as f:
        f.write("\n".join(lines) + "\n")
    print("\n".join(l for l in lines if l.startswith("#") or " 6 chip" in l))


if __name__ == "__main__":
    main()

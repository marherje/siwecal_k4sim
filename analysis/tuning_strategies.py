#!/usr/bin/env python
"""
sim_tuning 2: impact of each candidate digitiser change, emulated offline on the digitised simulation, against data.

Strategies (analysis/sim_tuning.py apply_*), each alone and combined:
  as now        final_v7 (or --campaign) as produced
  BL x1, x2     2.1 chip baseline drop B = slope x Q (Q = HG-ADC charge of the other hits of the chip): HG hits lose
                B, survive the per-slab trigger with P(a+B)/P(a). slope -0.95e-3 (cm_trigger_test e291+e287+e286 pooled), and twice
                that (upper bound: the measured B is biased up by sub-threshold deposits)
  BL chip       2.1 with the per-chip slope map (--slope-map npz: slope[slab, chip]) where measured
  LG(Q)         2.3 hits read in low gain scaled by the measured LG/HG plateau r(Q) / r(Q of the anchor fit)
                (--gain npz from tuning_gain.py, data of the group's set at 74 GeV)
  dropout       2.2 whole chips removed with the excess dropout probability of the TRAINING group (tuning_noise.py)
  cell eff      2.2 per-cell efficiency 2^-(local residual of log2 sim/data), capped at 1, from the TRAINING group
  all           BL x1 + LG(Q) + dropout + cell eff
Training group: --train (default th210_e20_26_-52); the others are validation (a map fitted on a sample always helps
that sample, so only the validation groups measure the strategy).

    python analysis/tuning_strategies.py <outdir> [--groups ...] [--train ...] [--campaign <dir>] [--tag v7]

Writes strategies.txt (figures of merit per group and strategy), layers_<group>.png (per-layer hits and E/hit
sim/data), chips_<group>.png (chip maps of log2 sim/data hits for every strategy), cells_<group>.png (local
residual cell maps, as now vs cell eff) and strategies.npz.
"""
import argparse
import os
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.colors import TwoSlopeNorm  # noqa: E402

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from analysis import sim_tuning as st  # noqa: E402

NL, NCHIP, NG = st.NL, st.NCHIP, st.NG
FOM = ["hits", "E", "eph", "layer_chi2", "chip_rms", "chip_eph_rms", "cell_rms", "cell_res_rms"]


def chip_layout(s):
    """(row, col) of every (slab, chip) on the 4x4 chip grid, from the mean cell of its hits."""
    rc = np.full((NL, NCHIP, 2), -1, int)
    for l in range(NL):
        for c in range(NCHIP):
            m = (s["slab"] == l) & (s["chip"] == c)
            if m.sum() > 20:
                rc[l, c] = (int(s["iy"][m].mean() // 8), int(s["ix"][m].mean() // 8))
    return rc


def chip_image(m, rc):
    img = np.full((NL, 4, 4), np.nan)
    for l in range(NL):
        for c in range(NCHIP):
            if rc[l, c, 0] >= 0:
                img[l, rc[l, c, 0], rc[l, c, 1]] = m[l, c]
    return img


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("outdir")
    ap.add_argument("--groups", default="th210_e20_26_-52,th210_e74_22_-49,th230_E4_e20,th230_E4_e74")
    ap.add_argument("--train", default="th210_e20_26_-52")
    ap.add_argument("--campaign", default=st.V7)
    ap.add_argument("--tag", default="v7")
    ap.add_argument("--gain", default=f"{st.OUT}/gain/gain_rq.npz")
    ap.add_argument("--noise", default=f"{st.OUT}/noise")
    ap.add_argument("--slope-map", default=f"{st.OUT}/occupancy/slope_map.npz")
    ap.add_argument("--slope", type=float, default=-0.95e-3)
    args = ap.parse_args()
    os.makedirs(args.outdir, exist_ok=True)
    rng = np.random.default_rng(7)
    groups = {g["name"]: g for g in st.load_groups()}
    names = args.groups.split(",")

    def samples(n):
        d = st.load(st.data_paths(groups[n]), cache=f"{st.OUT}/cache/data_{n}.npz")
        s = st.load(st.sim_paths(n, args.campaign), cache=f"{st.OUT}/cache/{args.tag}_{n}.npz")
        return d, s

    # ---- training-group maps
    dT, sT = samples(args.train)
    cT = st.compare(st.observables(sT), st.observables(dT))
    eff_cell = np.where(np.isfinite(cT["cell_res"]), np.minimum(2.0 ** (-cT["cell_res"]), 1.0), 1.0)
    nz = f"{args.noise}/noise_{args.train}.npz"
    p_drop = np.load(nz)["p_drop"] if os.path.exists(nz) else None
    slope_map = np.load(args.slope_map)["slope"] if os.path.exists(args.slope_map) else None
    gain = np.load(args.gain) if os.path.exists(args.gain) else None

    def lowgain(s, n):
        if gain is None:
            return None
        key = "th210_e74_22_-49_data" if "th210_e74_22_-49_data" in gain else [k for k in gain if k.endswith("_data")][0]
        q_c, r_q = gain["q_c"], gain[key]
        ok = np.isfinite(r_q)
        # the anchor (LG->HG line) was fitted on 99 GeV showers: reference = r at the median Q of LG-read hits there
        r_ref = float(np.interp(np.log(150.0), np.log(q_c[ok]), r_q[ok]))
        return st.apply_lowgain(s, q_c[ok], r_q[ok], r_ref)

    lines = [f"{'group':<20}{'strategy':<14}" + "".join(f"{k:>13}" for k in FOM)]
    store = {}
    for n in names:
        d, s = samples(n)
        th = groups[n]["set"]
        do = st.observables(d)
        strat = {"as now": s,
                 "BL x1": st.apply_baseline(s, th, args.slope, rng),
                 "BL x2": st.apply_baseline(s, th, 2 * args.slope, rng)}
        if slope_map is not None:
            strat["BL chip"] = st.apply_baseline(s, th, args.slope, rng, slope_map=slope_map)
        lg = lowgain(s, n)
        if lg is not None:
            strat["LG(Q)"] = lg
        if p_drop is not None:
            strat["dropout"] = st.apply_chip_dropout(s, p_drop, rng)
        strat["cell eff"] = st.apply_efficiency_map(s, eff_cell, rng)
        comb = st.apply_baseline(s, th, args.slope, rng)
        if lg is not None:
            comb = lowgain(comb, n)
        if p_drop is not None:
            comb = st.apply_chip_dropout(comb, p_drop, rng)
        strat["all"] = st.apply_efficiency_map(comb, eff_cell, rng)

        res = {k: st.compare(st.observables(v), do) for k, v in strat.items()}
        for k, c in res.items():
            lines.append(f"{n + (' [train]' if n == args.train else ''):<20}{k:<14}" +
                         "".join(f"{c[f]:13.3f}" for f in FOM))
            store[f"{n}|{k}"] = np.array([c[f] for f in FOM])
        print("\n".join(lines[-len(res):]), flush=True)

        # per layer
        fig, ax = plt.subplots(1, 2, figsize=(13, 4.3))
        for i, (k, c) in enumerate(res.items()):
            ls = "-" if k == "as now" else "--"
            ax[0].plot(np.arange(NL), c["layer_ratio"], ls, marker="o", ms=3, label=k, color=plt.cm.tab10(i))
            ax[1].plot(np.arange(NL), c["layer_eph_ratio"], ls, marker="o", ms=3, color=plt.cm.tab10(i))
        for a, t in zip(ax, ("hits / event", "E / hit")):
            a.axhline(1, color="k", lw=0.8); a.set_xlabel("layer"); a.set_ylabel(f"{t}, sim / data")
        ax[0].legend(fontsize=8, ncol=2); fig.suptitle(n); fig.tight_layout()
        fig.savefig(f"{args.outdir}/layers_{n}.png", dpi=110); plt.close(fig)

        # chip maps: one row per strategy, 15 slabs
        rc = chip_layout(d)
        ks = list(res)
        fig, ax = plt.subplots(len(ks), NL, figsize=(1.25 * NL, 1.35 * len(ks) + 0.6), squeeze=False)
        norm = TwoSlopeNorm(vmin=-0.8, vcenter=0, vmax=0.8)
        for i, k in enumerate(ks):
            img = chip_image(res[k]["chip_log2"], rc)
            for l in range(NL):
                im = ax[i, l].imshow(img[l], origin="lower", cmap="RdBu_r", norm=norm)
                ax[i, l].set_xticks([]); ax[i, l].set_yticks([])
                if i == 0:
                    ax[i, l].set_title(f"L{l}", fontsize=8)
            ax[i, 0].set_ylabel(f"{k}\nrms {res[k]['chip_rms']:.3f}", fontsize=8)
        fig.colorbar(im, ax=ax, shrink=0.6, label="log2 hits sim/data")
        fig.suptitle(f"{n}: per-chip hits/event, sim / data")
        fig.savefig(f"{args.outdir}/chips_{n}.png", dpi=100, bbox_inches="tight"); plt.close(fig)

        # cell maps: local residual, as now vs cell eff
        fig, ax = plt.subplots(2, NL, figsize=(1.25 * NL, 3.4), squeeze=False)
        for i, k in enumerate(("as now", "cell eff")):
            for l in range(NL):
                im = ax[i, l].imshow(res[k]["cell_res"][l], origin="lower", cmap="RdBu_r", vmin=-1, vmax=1)
                ax[i, l].set_xticks([]); ax[i, l].set_yticks([])
                if i == 0:
                    ax[i, l].set_title(f"L{l}", fontsize=8)
            ax[i, 0].set_ylabel(f"{k}\nrms {res[k]['cell_res_rms']:.3f}", fontsize=8)
        fig.colorbar(im, ax=ax, shrink=0.7, label="log2 sim/data - 5x5 median")
        fig.suptitle(f"{n}: local residual per cell")
        fig.savefig(f"{args.outdir}/cells_{n}.png", dpi=100, bbox_inches="tight"); plt.close(fig)

    open(f"{args.outdir}/strategies.txt", "w").write("\n".join(lines) + "\n")
    np.savez(f"{args.outdir}/strategies.npz", fom=np.array(FOM), **store)


if __name__ == "__main__":
    main()

#!/usr/bin/env python
"""
sim_tuning 2.1b: fit the chip-occupancy effect from the hits themselves, and look at where nhits and E/hit differ.

The baseline emulation (sim_tuning.apply_baseline: B = k x max(Q - Q0, 0), Q = HG-ADC charge of the other hits of
the chip; HG hits lose B and survive the per-slab turn-on with P(a+B)/P(a)) is scanned over
  lin    k in KGRID, Q0 = 0                         (k = -0.95e-3 is the pedestal measurement of cm_trigger_test)
  onset  k in KGRID, Q0 in Q0GRID                   (no shift below a chip charge Q0)
  map    the measured per-chip slope map x f, f in FGRID
on the TRAINING group (--train, default the highest energy, where the occupancy is highest), with the distribution of
the number of hits per chip and event (chip multiplicity) as the figure of merit:
  mult_chi2   sum over multiplicity bins of (n_sim - n_data)^2 / (n_sim/N_sim + n_data/N_data), per bin, n = chip-events
              per event in the bin;
and validated on the other groups. Also per model: hits/event, E/hit (sim/data), layer chi2 of hits.

Diagnostics (as now / measured map / best fit): chip multiplicity sim/data, hit-energy spectrum sim/data,
E/hit sim/data in chip-multiplicity bins and in radial bins about the event barycentre (core / mid / halo).

    python analysis/tuning_occupancy.py <outdir> [--campaign <dir> --tag h2line] [--groups ...] [--train ...]

Writes occupancy_fit.txt, occupancy_scan.png, occupancy_mult.png, occupancy_espec.png, occupancy_eph.png,
occupancy_fit.npz.
"""
import argparse
import os
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from analysis import sim_tuning as st  # noqa: E402

KGRID = np.array([0.0, -0.5, -0.95, -1.5, -2.0, -2.5, -3.0, -4.0]) * 1e-3
Q0GRID = np.array([2000.0, 4000.0, 6000.0])
FGRID = np.array([1.0, 1.5, 2.0, 2.5, 3.0])
MB = np.array([1, 6, 11, 16, 21, 26, 31, 36, 41, 46, 65])          # chip multiplicity bins
EB = np.array([0.5, 1.2, 2, 3, 5, 10, 20, 50, 1e9])                 # hit energy [MIP]
RB = np.array([0, 11, 33, 1e9])                                     # radius about the barycentre [mm]
RLAB = ["core r<11", "mid 11-33", "halo r>33"]


def mult(s):
    _, c = np.unique(st.chip_key(s), return_counts=True)
    return np.histogram(c, MB)[0].astype(float)


def mult_chi2(s, d):
    ns, nd = s["nev"], d["nev"]
    a, b = mult(s), mult(d)
    v = a / ns ** 2 + b / nd ** 2
    ok = (a + b) > 20
    return float(np.mean((a[ok] / ns - b[ok] / nd) ** 2 / v[ok]))


def eph_bins(s):
    occ = st.chip_nhits(s)
    r = np.hypot(s["x"] - s["bar_x"][s["ev"]], s["y"] - s["bar_y"][s["ev"]])
    e = s["e"].astype(float)
    io = np.clip(np.digitize(occ, MB) - 1, 0, len(MB) - 2)
    ir = np.clip(np.digitize(r, RB) - 1, 0, len(RB) - 2)
    ho = np.bincount(io, minlength=len(MB) - 1); eo = np.bincount(io, weights=e, minlength=len(MB) - 1)
    hr = np.bincount(ir, minlength=len(RB) - 1); er = np.bincount(ir, weights=e, minlength=len(RB) - 1)
    return {"occ_eph": eo / np.maximum(ho, 1), "occ_n": ho / s["nev"], "r_eph": er / np.maximum(hr, 1),
            "r_n": hr / s["nev"], "espec": np.histogram(e, EB)[0] / s["nev"]}


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("outdir")
    ap.add_argument("--campaign", default=f"{st.OUT}/upstream/h2line")
    ap.add_argument("--tag", default="h2line")
    ap.add_argument("--groups", default="th230_E4_e10,th230_E4_e20,th230_E4_e52,th230_E4_e74,th210_e20_26_-52")
    ap.add_argument("--train", default="th230_E4_e74")
    ap.add_argument("--slope-map", default=f"{st.OUT}/occupancy/slope_map.npz")
    args = ap.parse_args()
    os.makedirs(args.outdir, exist_ok=True)
    groups = {g["name"]: g for g in st.load_groups()}
    names = args.groups.split(",")
    smap = np.load(args.slope_map)["slope"]

    models = [("lin", k, 0.0, None) for k in KGRID]
    models += [("onset", k, q0, None) for q0 in Q0GRID for k in KGRID[1:]]
    models += [("map", -0.95e-3, 0.0, f) for f in FGRID]
    mname = lambda m: (f"lin k={m[1] * 1e3:.2f}" if m[0] == "lin" else
                       f"onset k={m[1] * 1e3:.2f} Q0={m[2]:.0f}" if m[0] == "onset" else f"map x{m[3]:.1f}")

    def apply(s, th, m, seed=11):
        rng = np.random.default_rng(seed)
        if m[0] == "lin" and m[1] == 0.0:
            return s
        if m[0] == "map":
            return st.apply_baseline(s, th, m[1], rng, slope_map=smap * m[3])
        return st.apply_baseline(s, th, m[1], rng, q0_adc=m[2])

    data, sims = {}, {}
    for n in names:
        data[n] = st.load(st.data_paths(groups[n]), cache=f"{st.OUT}/cache/data_{n}.npz")
        sims[n] = st.load(st.sim_paths(n, args.campaign), cache=f"{st.OUT}/cache/{args.tag}_{n}.npz")
    dobs = {n: st.observables(data[n]) for n in names}

    # ---- scan every model on every group
    R = {}
    for m in models:
        for n in names:
            s = apply(sims[n], groups[n]["set"], m)
            c = st.compare(st.observables(s), dobs[n])
            R[(mname(m), n)] = {"mult_chi2": mult_chi2(s, data[n]), "hits": c["hits"], "eph": c["eph"],
                                "E": c["E"], "layer_chi2": c["layer_chi2"]}
        print(mname(m), " ".join(f"{n}:{R[(mname(m), n)]['mult_chi2']:.0f}/{R[(mname(m), n)]['hits']:.3f}"
                                 for n in names), flush=True)
    best = min(models, key=lambda m: R[(mname(m), args.train)]["mult_chi2"])
    bestlin = min([m for m in models if m[0] == "lin"], key=lambda m: R[(mname(m), args.train)]["mult_chi2"])
    bestmap = min([m for m in models if m[0] == "map"], key=lambda m: R[(mname(m), args.train)]["mult_chi2"])
    ref = {"as now": models[0], "measured map": ("map", -0.95e-3, 0.0, 1.0), "best lin": bestlin,
           "best map": bestmap, "best": best}

    lines = [f"training group {args.train}; figure of merit: chip multiplicity chi2 per bin",
             f"best overall: {mname(best)}; best linear: {mname(bestlin)}; best scaled map: {mname(bestmap)}", "",
             f"{'model':<26}" + "".join(f"{n:>34}" for n in names),
             f"{'':<26}" + "".join(f"{'multchi2 hits E/hit layerchi2':>34}" for _ in names)]
    for m in models:
        lines.append(f"{mname(m):<26}" + "".join(
            f"{R[(mname(m), n)]['mult_chi2']:9.0f} {R[(mname(m), n)]['hits']:6.3f} {R[(mname(m), n)]['eph']:6.3f}"
            f" {R[(mname(m), n)]['layer_chi2']:9.0f}" for n in names))
    open(f"{args.outdir}/occupancy_fit.txt", "w").write("\n".join(lines) + "\n")
    print("\n".join(lines[:2]))

    # ---- scan figure (linear and onset models; map as a separate x axis)
    fig, ax = plt.subplots(1, 3, figsize=(16, 4.3))
    for j, n in enumerate(names):
        col = plt.cm.viridis(j / max(len(names) - 1, 1))
        lw = 2.2 if n == args.train else 1.2
        ks = -KGRID * 1e3
        ax[0].plot(ks, [R[(mname(("lin", k, 0.0, None)), n)]["mult_chi2"] for k in KGRID], "-o", ms=3, color=col,
                   lw=lw, label=n + (" [train]" if n == args.train else ""))
        ax[1].plot(ks, [R[(mname(("lin", k, 0.0, None)), n)]["hits"] for k in KGRID], "-o", ms=3, color=col, lw=lw)
        ax[2].plot(ks, [R[(mname(("lin", k, 0.0, None)), n)]["eph"] for k in KGRID], "-o", ms=3, color=col, lw=lw)
    for a in ax:
        a.axvline(0.95, color="k", ls=":", lw=1)
        a.set_xlabel("-k  [ADC baseline shift per 1000 ADC of chip charge]")
    ax[0].set_yscale("log"); ax[0].set_ylabel("chip multiplicity chi2 / bin"); ax[0].legend(fontsize=8)
    ax[1].axhline(1, color="k", lw=0.8); ax[1].set_ylabel("hits / event, sim / data")
    ax[2].axhline(1, color="k", lw=0.8); ax[2].set_ylabel("E / hit, sim / data")
    fig.suptitle("Linear chip-baseline model B = k Q: scan (dotted: pedestal measurement)")
    fig.tight_layout(); fig.savefig(f"{args.outdir}/occupancy_scan.png", dpi=110); plt.close(fig)

    # ---- diagnostics for the reference models
    D = {}
    for lab, m in ref.items():
        for n in names:
            s = apply(sims[n], groups[n]["set"], m)
            D[(lab, n)] = {"mult": mult(s) / s["nev"], **eph_bins(s)}
    for n in names:
        D[("data", n)] = {"mult": mult(data[n]) / data[n]["nev"], **eph_bins(data[n])}
    labs = ["as now", "measured map", "best"]
    cols = {"as now": "#8c8c8c", "measured map": "#eda100", "best": "#c0392b"}
    mc = 0.5 * (MB[1:] + MB[:-1]) - 0.5

    fig, ax = plt.subplots(1, len(names), figsize=(4.0 * len(names), 3.8), squeeze=False)
    for j, n in enumerate(names):
        for lab in labs:
            with np.errstate(divide="ignore", invalid="ignore"):
                ax[0, j].plot(mc, D[(lab, n)]["mult"] / D[("data", n)]["mult"], "-o", ms=3, color=cols[lab],
                              label=lab if lab != "best" else f"best: {mname(best)}")
        ax[0, j].axhline(1, color="k", lw=0.8); ax[0, j].set_ylim(0.4, 2.5); ax[0, j].set_yscale("log")
        ax[0, j].set_title(n, fontsize=10); ax[0, j].set_xlabel("hits in the chip (event)")
    ax[0, 0].set_ylabel("chip-events, sim / data"); ax[0, 0].legend(fontsize=7)
    fig.tight_layout(); fig.savefig(f"{args.outdir}/occupancy_mult.png", dpi=110); plt.close(fig)

    ec = np.sqrt(EB[:-1] * np.minimum(EB[1:], 200))
    fig, ax = plt.subplots(1, len(names), figsize=(4.0 * len(names), 3.8), squeeze=False)
    for j, n in enumerate(names):
        for lab in labs:
            ax[0, j].plot(ec, D[(lab, n)]["espec"] / np.maximum(D[("data", n)]["espec"], 1e-12), "-o", ms=3,
                          color=cols[lab], label=lab)
        ax[0, j].axhline(1, color="k", lw=0.8); ax[0, j].set_xscale("log"); ax[0, j].set_ylim(0.6, 1.6)
        ax[0, j].set_title(n, fontsize=10); ax[0, j].set_xlabel("hit energy [MIP]")
    ax[0, 0].set_ylabel("hits per energy bin, sim / data"); ax[0, 0].legend(fontsize=7)
    fig.tight_layout(); fig.savefig(f"{args.outdir}/occupancy_espec.png", dpi=110); plt.close(fig)

    fig, ax = plt.subplots(2, len(names), figsize=(4.0 * len(names), 7.0), squeeze=False)
    for j, n in enumerate(names):
        dd = D[("data", n)]
        for lab in labs:
            x = D[(lab, n)]
            ok = dd["occ_n"] * data[n]["nev"] > 200
            ax[0, j].plot(mc[ok], (x["occ_eph"] / dd["occ_eph"])[ok], "-o", ms=3, color=cols[lab], label=lab)
            ax[1, j].plot(np.arange(3), x["r_eph"] / dd["r_eph"], "-o", ms=4, color=cols[lab])
            for k in range(3):
                ax[1, j].annotate(f"{x['r_n'][k] / dd['r_n'][k]:.3f}", (k, x["r_eph"][k] / dd["r_eph"][k]),
                                  fontsize=6, color=cols[lab], xytext=(3, 3), textcoords="offset points")
        for a in ax[:, j]:
            a.axhline(1, color="k", lw=0.8)
        ax[0, j].set_title(n, fontsize=10); ax[0, j].set_xlabel("hits in the chip (event)")
        ax[1, j].set_xticks(range(3)); ax[1, j].set_xticklabels(RLAB, fontsize=8)
    ax[0, 0].set_ylabel("E / hit, sim / data"); ax[1, 0].set_ylabel("E / hit, sim / data (labels: hits sim/data)")
    ax[0, 0].legend(fontsize=7)
    fig.tight_layout(); fig.savefig(f"{args.outdir}/occupancy_eph.png", dpi=110); plt.close(fig)

    np.savez(f"{args.outdir}/occupancy_fit.npz", best=np.array([best[0], best[1], best[2], best[3] or 0], dtype=object),
             **{f"{k[0]}|{k[1]}": np.array([v["mult_chi2"], v["hits"], v["eph"], v["E"], v["layer_chi2"]])
                for k, v in R.items()})
    print(f"\nreference models: { {k: mname(v) for k, v in ref.items()} }")
    for n in names:
        print(n, " ".join(f"{lab}: hits {R[(mname(m), n)]['hits']:.3f} eph {R[(mname(m), n)]['eph']:.3f} "
                          f"mult {R[(mname(m), n)]['mult_chi2']:.0f} |" for lab, m in ref.items()))
        for lab in ["as now", "best"]:
            x, dd = D[(lab, n)], D[("data", n)]
            print(f"   {lab}: E/hit by radius core/mid/halo " + " ".join(f"{a:.3f}" for a in x["r_eph"] / dd["r_eph"])
                  + "  hits " + " ".join(f"{a:.3f}" for a in x["r_n"] / dd["r_n"]))


if __name__ == "__main__":
    main()

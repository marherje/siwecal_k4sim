#!/usr/bin/env python
"""
Longitudinal shower shape, data against one or more simulations of the same beam points: does the simulated shower
start (and peak) where the data's does? The figure of merit for the amount of material upstream of the calorimeter.

    python analysis/upstream_profile.py <outdir> --sim "v7=<campaign dir>" [--sim "Al5=<dir>" ...] \
        [--groups th210_e20_26_-52,th230_E4_e20,...] [--kind digi]

<campaign dir> holds pid/<group>_<kind>/*.valtree.root (make_sim_dag.py layout). Data: Reconstruction valtrees of the
group's runs (analysis/beam_points.txt). Shower-like events (nhit > 0.5 x p90, as final_panels.py), 0.5 MIP hit cut
of the PID stage in both.

Per group and sample:
  f_l      mean fraction of the event energy in layer l (normalised profile: insensitive to the overall scale)
  zbary    mean energy barycentre [layer]
  front    fraction of the energy in layers 0-2
  chi2     sum over layers of (f_sim - f_data)^2 / (var_sim + var_data), per degree of freedom (15 layers), with
           var_sim rescaled to the data's number of events
Writes <outdir>/upstream_profile.txt, profile_<group>.png (normalised profiles + sim/data per layer for every
simulation) and upstream_summary.png (zbary, front fraction and chi2 against the beam energy, one line per sim),
and a cache per (sample, group) in <outdir>/cache/.
"""
import argparse
import glob
import os
import sys

import numpy as np
import uproot
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(_HERE))
from analysis.final_panels import load_groups, RECO  # noqa: E402

NL = 15
COLORS = ["#c0392b", "#eda100", "#1baf7a", "#4a3aa7", "#8c8c8c", "#e377c2", "#17becf", "#7f4f24"]


def sample_profile(paths, cache, nmax=40000):
    if os.path.exists(cache):
        return dict(np.load(cache))
    NH, EPL, Z = [], [], []
    left = nmax
    for p in paths:
        if left <= 0:
            break
        a = uproot.open(p)["ecal"].arrays(["nhit_chan", "energy_per_layer", "zbary"], entry_stop=left, library="np")
        left -= len(a["nhit_chan"])
        NH.append(a["nhit_chan"]); EPL.append(np.stack(a["energy_per_layer"]).astype(float)); Z.append(a["zbary"])
    nh, epl, zb = np.concatenate(NH), np.concatenate(EPL), np.concatenate(Z)
    k = nh > 0.5 * np.percentile(nh, 90)
    epl, zb = epl[k], zb[k]
    tot = epl.sum(1)
    f = epl / tot[:, None]
    out = {"f": f.mean(0), "fvar": f.var(0) / len(f), "eabs": epl.mean(0), "zbary": np.nanmean(zb),
           "zbary_err": np.nanstd(zb) / np.sqrt(len(zb)), "front": f[:, :3].sum(1).mean(), "n": len(f)}
    os.makedirs(os.path.dirname(cache), exist_ok=True)
    np.savez(cache, **out)
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("outdir")
    ap.add_argument("--sim", action="append", required=True, help="label=<campaign dir>")
    ap.add_argument("--groups", default="th210_e7.5_17_-50,th210_e10_19_-51,th210_e20_26_-52,th210_e34_27_-53,"
                    "th210_e52_24_-51,th210_e74_22_-49,th230_E4_e7.5,th230_E4_e10,th230_E4_e20,th230_E4_e34,"
                    "th230_E4_e52,th230_E4_e74")
    ap.add_argument("--kind", default="digi")
    ap.add_argument("--max-events", type=int, default=40000)
    args = ap.parse_args()
    os.makedirs(args.outdir, exist_ok=True)
    sims = [s.split("=", 1) for s in args.sim]
    groups = {g["name"]: g for g in load_groups()}
    names = [n for n in args.groups.split(",") if n in groups]

    res = {}
    lines = ["group sample n zbary zbary_err front chi2/ndf  f_0..f_14"]
    for name in names:
        g = groups[name]
        dpaths = [p for r in g["runs"] for p in glob.glob(f"{RECO}/{r}/*.valtree.root")]
        if not dpaths:
            print("no data for", name); continue
        res[name] = {"data": sample_profile(dpaths, f"{args.outdir}/cache/data/{name}.npz", args.max_events)}
        for lab, d in sims:
            sp = sorted(glob.glob(f"{d}/pid/{name}_{args.kind}/*.valtree.root"))
            if sp:
                res[name][lab] = sample_profile(sp, f"{args.outdir}/cache/{lab}/{name}_{args.kind}.npz", args.max_events)
        dd = res[name]["data"]
        for lab, r in res[name].items():
            # sim variance taken at the data's statistics, so samples of different size compare on the same footing
            r["chi2"] = 0.0 if lab == "data" else float(
                ((r["f"] - dd["f"]) ** 2 / (r["fvar"] * r["n"] / dd["n"] + dd["fvar"])).sum() / NL)
            lines.append(f"{name} {lab} {int(r['n'])} {r['zbary']:.3f} {r['zbary_err']:.3f} {r['front']:.4f} "
                         f"{r['chi2']:.1f}  " + " ".join(f"{x:.4f}" for x in r["f"]))

        fig, ax = plt.subplots(2, 1, figsize=(7, 6.5), sharex=True, gridspec_kw={"height_ratios": [2, 1]})
        L = np.arange(NL)
        ax[0].errorbar(L, dd["f"], np.sqrt(dd["fvar"]), fmt="o", color="#29a3dc", label="Data")
        for i, (lab, _) in enumerate(sims):
            if lab not in res[name]:
                continue
            r = res[name][lab]
            ax[0].plot(L, r["f"], "-", color=COLORS[i % len(COLORS)], label=f"{lab} (chi2/ndf {r['chi2']:.0f})")
            ax[1].plot(L, r["f"] / dd["f"], "-o", ms=3, color=COLORS[i % len(COLORS)])
        ax[1].axhline(1, color="k", lw=0.8)
        ax[0].set_ylabel("mean energy fraction per layer"); ax[1].set_ylabel("sim / data"); ax[1].set_xlabel("layer")
        ax[0].set_title(name); ax[0].legend(fontsize=8)
        fig.tight_layout(); fig.savefig(f"{args.outdir}/profile_{name}.png", dpi=110); plt.close(fig)

    open(f"{args.outdir}/upstream_profile.txt", "w").write("\n".join(lines) + "\n")
    print("\n".join(l[:80] for l in lines))

    # summary against energy, one panel row per threshold set
    sets = sorted({groups[n]["set"] for n in res})
    fig, ax = plt.subplots(len(sets), 3, figsize=(14, 3.6 * len(sets)), squeeze=False)
    for si, s in enumerate(sets):
        ns = [n for n in res if groups[n]["set"] == s]
        E = np.array([float(groups[n]["E"]) for n in ns])
        o = np.argsort(E)
        for i, (lab, _) in enumerate(sims):
            ok = [n for n in np.array(ns)[o] if lab in res[n]]
            if not ok:
                continue
            Ek = [float(groups[n]["E"]) for n in ok]
            c = COLORS[i % len(COLORS)]
            ax[si, 0].plot(Ek, [res[n][lab]["zbary"] - res[n]["data"]["zbary"] for n in ok], "-o", color=c, label=lab)
            ax[si, 1].plot(Ek, [res[n][lab]["front"] / res[n]["data"]["front"] for n in ok], "-o", color=c)
            ax[si, 2].plot(Ek, [res[n][lab]["chi2"] for n in ok], "-o", color=c)
        ax[si, 0].axhline(0, color="k", lw=0.8); ax[si, 1].axhline(1, color="k", lw=0.8)
        ax[si, 0].set_ylabel(f"{s}: zbary sim - data [layer]"); ax[si, 1].set_ylabel("front (L0-2) sim/data")
        ax[si, 2].set_ylabel("profile chi2 / ndf"); ax[si, 2].set_yscale("log")
        for a in ax[si]:
            a.set_xlabel("beam energy [GeV]")
        ax[si, 0].legend(fontsize=8)
    fig.tight_layout(); fig.savefig(f"{args.outdir}/upstream_summary.png", dpi=110); plt.close(fig)


if __name__ == "__main__":
    main()

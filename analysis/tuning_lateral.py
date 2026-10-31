#!/usr/bin/env python
"""
sim_tuning 2.4: is the data shower laterally wider (more leakage out of the core) than the simulation?

Per group, data against one or more simulations (label=campaign dir, digitised; 'nodigi' kind with label:nodigi):
  radial      energy fraction against the distance r to the event barycentre, in layer blocks 0-2, 3-5, 6-8, 9-14,
              and sim/data
  R90, f>22   per layer, the radius holding 90 % of the layer's mean energy, and the fraction beyond 22 mm (4 pads)
  wander      per layer, median distance between the layer's own energy barycentre and the event barycentre
  map         energy density on the pad grid about the event barycentre (+-15 pads), layers 3-8, log2 sim/data

    python analysis/tuning_lateral.py <outdir> --sim "v7=<dir>" [--sim "v7 nodigi=<dir>:nodigi"] [--groups ...]

Writes lateral_<group>.png, lateral.txt.
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

NL = st.NL
BLOCKS = [(0, 2), (3, 5), (6, 8), (9, 14)]
RB = np.array([0, 2.8, 5.5, 8.3, 11, 16.5, 22, 27.5, 33, 44, 55, 66, 88, 110, 150.0])
H = 15


def lateral(s):
    nev = int(s["nev"])
    ev, e, sl = s["ev"], s["e"].astype(float), s["slab"].astype(int)
    dx, dy = s["x"] - s["bar_x"][ev], s["y"] - s["bar_y"][ev]
    r = np.hypot(dx, dy)
    rb = np.clip(np.digitize(r, RB) - 1, 0, len(RB) - 2)
    prof = np.bincount(sl * (len(RB) - 1) + rb, weights=e, minlength=NL * (len(RB) - 1)).reshape(NL, -1) / nev
    # layer barycentre wander
    key = ev.astype(np.int64) * NL + sl
    u, inv = np.unique(key, return_inverse=True)
    w = np.bincount(inv, weights=e)
    lx = np.bincount(inv, weights=e * s["x"]) / w; ly = np.bincount(inv, weights=e * s["y"]) / w
    uev, ul = u // NL, u % NL
    wd = np.hypot(lx - s["bar_x"][uev], ly - s["bar_y"][uev])
    wander = np.array([np.median(wd[(ul == l) & (w > 3)]) if ((ul == l) & (w > 3)).sum() > 50 else np.nan
                       for l in range(NL)])
    # pad-grid density about the barycentre, layers 3-8
    m = (sl >= 3) & (sl <= 8)
    ix = np.rint(dx[m] / st.PITCH).astype(int) + H; iy = np.rint(dy[m] / st.PITCH).astype(int) + H
    ok = (ix >= 0) & (ix <= 2 * H) & (iy >= 0) & (iy <= 2 * H)
    dens = np.zeros((2 * H + 1, 2 * H + 1))
    np.add.at(dens, (iy[ok], ix[ok]), e[m][ok])
    return {"prof": prof, "wander": wander, "dens": dens / nev}


def r90(p):
    c = np.cumsum(p) / max(p.sum(), 1e-12)
    return float(np.interp(0.9, c, RB[1:]))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("outdir")
    ap.add_argument("--sim", action="append", required=True)
    ap.add_argument("--groups", default="th210_e20_26_-52,th210_e74_22_-49,th230_E4_e20,th230_E4_e74")
    args = ap.parse_args()
    os.makedirs(args.outdir, exist_ok=True)
    groups = {g["name"]: g for g in st.load_groups()}
    lines = ["group sample  R90[mm] per layer | f(r>22mm) per layer | wander[mm] per layer"]
    for n in args.groups.split(","):
        g = groups[n]
        d = st.load(st.data_paths(g), cache=f"{st.OUT}/cache/data_{n}.npz")
        samples = {"data": lateral(d)}
        for spec in args.sim:
            lab, path = spec.split("=", 1)
            kind = "digi"
            if path.endswith(":nodigi"):
                path, kind = path[:-7], "nodigi"
            tag = lab.replace(" ", "_")
            paths = st.sim_paths(n, path, kind)
            if paths:
                samples[lab] = lateral(st.load(paths, cache=f"{st.OUT}/cache/{tag}_{n}.npz"))
        for lab, L in samples.items():
            R = [r90(L["prof"][l]) for l in range(NL)]
            f22 = [L["prof"][l][RB[1:] > 22].sum() / max(L["prof"][l].sum(), 1e-12) for l in range(NL)]
            lines.append(f"{n} {lab} " + " ".join(f"{x:.1f}" for x in R) + " | " + " ".join(f"{x:.3f}" for x in f22)
                         + " | " + " ".join(f"{x:.1f}" for x in L["wander"]))
        print("\n".join(lines[-len(samples):]), flush=True)

        D = samples["data"]
        fig, ax = plt.subplots(2, 4, figsize=(17, 7.5))
        rc = 0.5 * (RB[1:] + RB[:-1])
        for b, (l0, l1) in enumerate(BLOCKS):
            pd_ = D["prof"][l0:l1 + 1].sum(0)
            ax[0, b].step(rc, pd_ / pd_.sum(), where="mid", color="k", label="data")
            for i, (lab, L) in enumerate(samples.items()):
                if lab == "data":
                    continue
                ps = L["prof"][l0:l1 + 1].sum(0)
                ax[0, b].step(rc, ps / ps.sum(), where="mid", color=plt.cm.tab10(i), label=lab)
                ax[1, b].plot(rc, (ps / ps.sum()) / np.maximum(pd_ / pd_.sum(), 1e-12), "-o", ms=3,
                              color=plt.cm.tab10(i))
            ax[0, b].set_yscale("log"); ax[0, b].set_title(f"layers {l0}-{l1}")
            ax[1, b].axhline(1, color="k", lw=0.8); ax[1, b].set_xscale("log"); ax[0, b].set_xscale("log")
            ax[1, b].set_xlabel("r to the event barycentre [mm]"); ax[1, b].set_ylim(0.4, 1.6)
        ax[0, 0].set_ylabel("energy fraction"); ax[1, 0].set_ylabel("sim / data"); ax[0, 0].legend(fontsize=8)
        fig.suptitle(f"{n}: lateral energy profile"); fig.tight_layout()
        fig.savefig(f"{args.outdir}/lateral_{n}.png", dpi=110); plt.close(fig)

        sims = [k for k in samples if k != "data"]
        fig, ax = plt.subplots(1, len(sims) + 1, figsize=(4.4 * (len(sims) + 1), 4))
        ax = np.atleast_1d(ax)
        for i, lab in enumerate(sims):
            with np.errstate(divide="ignore", invalid="ignore"):
                rr = np.log2(samples[lab]["dens"] / D["dens"])
            im = ax[i].imshow(rr, origin="lower", cmap="RdBu_r", vmin=-1, vmax=1,
                              extent=[-H - 0.5, H + 0.5, -H - 0.5, H + 0.5])
            ax[i].set_title(f"{lab}: log2 E density sim/data, L3-8", fontsize=9)
            ax[i].set_xlabel("dx [pads]")
            ax[-1].plot(np.arange(NL), samples[lab]["wander"], "-o", ms=3, label=lab, color=plt.cm.tab10(i + 1))
        fig.colorbar(im, ax=ax[:-1])
        ax[-1].plot(np.arange(NL), D["wander"], "-o", color="k", label="data")
        ax[-1].set_xlabel("layer"); ax[-1].set_ylabel("median |layer - event barycentre| [mm]"); ax[-1].legend(fontsize=8)
        fig.suptitle(n)
        fig.savefig(f"{args.outdir}/lateralmap_{n}.png", dpi=110, bbox_inches="tight"); plt.close(fig)
    open(f"{args.outdir}/lateral.txt", "w").write("\n".join(lines) + "\n")


if __name__ == "__main__":
    main()

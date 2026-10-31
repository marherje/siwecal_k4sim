#!/usr/bin/env python
"""
sim_tuning 1: best amount of extra upstream material from the Al scan of upstream_profile.py.

    python analysis/upstream_fit.py <scan dir> [--groups ...]

Reads <scan dir>/cache/{data,v7,up0pXX,wd10X}/<group>_digi.npz. Per group, against the extra Al thickness t [X0]:
  chi2(t)    parabola through the 3 lowest points -> t_min, and the range where chi2 < 1.1 chi2_min (the profile
             chi2 is dominated by the remaining shape mismatch, not by statistics, so Delta chi2 = 1 is meaningless)
  zbary(t)   linear fit -> t where zbary_sim = zbary_data
  front(t)   linear fit -> t where front_sim = front_data
Writes upstream_fit.txt and upstream_fit.png (chi2, zbary and front against t per group, the W-density variants as
open markers at the t of equal zbary shift; per-layer sim/data at the best t).
"""
import argparse
import os
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

AL = {"v7": 0.0, "up0p05": 0.05, "up0p10": 0.10, "up0p15": 0.15, "up0p20": 0.20, "up0p30": 0.30}
WD = {"wd103": "W x1.03", "wd105": "W x1.05"}
NL = 15


def chi2(r, d):
    return float(((r["f"] - d["f"]) ** 2 / (r["fvar"] * r["n"] / d["n"] + d["fvar"])).sum() / NL)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("scan")
    ap.add_argument("--groups", default="th210_e20_26_-52,th230_E4_e10,th230_E4_e20,th230_E4_e52,th230_E4_e74")
    args = ap.parse_args()
    names = args.groups.split(",")
    ld = lambda lab, n: dict(np.load(f"{args.scan}/cache/{lab}/{n}{'' if lab == 'data' else '_digi'}.npz"))
    lines = ["group  t_chi2min  t(chi2<1.1min)  t_zbary  t_front  chi2(v7)  chi2(tmin)  chi2(W1.03)  chi2(W1.05)"]
    fig, ax = plt.subplots(1, 4, figsize=(19, 4.2))
    best = {}
    for i, n in enumerate(names):
        d = ld("data", n)
        R = {k: ld(k, n) for k in list(AL) + list(WD)}
        t = np.array(list(AL.values()))
        c = np.array([chi2(R[k], d) for k in AL])
        z = np.array([R[k]["zbary"] - d["zbary"] for k in AL])
        f = np.array([R[k]["front"] / d["front"] for k in AL])
        j = np.argsort(c)[:3]
        p = np.polyfit(t[j], c[j], 2)
        tmin = -p[1] / (2 * p[0]) if p[0] > 0 else t[j[0]]
        cmin = np.polyval(p, tmin)
        tt = np.linspace(0, 0.3, 301)
        ok = tt[np.polyval(np.polyfit(t, c, 3), tt) < 1.1 * cmin]
        tz = np.roots(np.polyfit(t, z, 1))[0]
        tf = np.roots(np.polyfit(t, f - 1, 1))[0]
        best[n] = tmin
        lines.append(f"{n} {tmin:.3f} {ok.min() if len(ok) else np.nan:.2f}-{ok.max() if len(ok) else np.nan:.2f} "
                     f"{tz:.3f} {tf:.3f} {c[0]:.0f} {cmin:.0f} {chi2(R['wd103'], d):.0f} {chi2(R['wd105'], d):.0f}")
        col = plt.cm.viridis(i / max(len(names) - 1, 1))
        ax[0].plot(t, c / c[0], "o-", color=col, label=n)
        ax[0].plot(tt, np.polyval(p, tt) / c[0], ":", color=col)
        ax[1].plot(t, z, "o-", color=col); ax[2].plot(t, f, "o-", color=col)
        for k, m in (("wd103", "s"), ("wd105", "D")):
            # place the W variant at the Al thickness giving the same zbary shift
            tw = np.interp(-(R[k]["zbary"] - d["zbary"] - z[0]), -(z - z[0]), t)
            ax[0].plot(tw, chi2(R[k], d) / c[0], m, mfc="none", color=col, ms=7)
            ax[2].plot(tw, R[k]["front"] / d["front"], m, mfc="none", color=col, ms=7)
        kb = min(AL, key=lambda k: abs(AL[k] - tmin))
        ax[3].plot(np.arange(NL), R["v7"]["f"] / d["f"], "--", color=col, lw=1)
        ax[3].plot(np.arange(NL), R[kb]["f"] / d["f"], "-o", ms=3, color=col, label=f"{n}: {kb}")
    ax[0].set_ylabel("profile chi2 / chi2(v7)"); ax[0].legend(fontsize=7)
    ax[0].plot([], [], "s", mfc="none", color="k", label="W x1.03 (at equal zbary shift)")
    ax[0].plot([], [], "D", mfc="none", color="k", label="W x1.05"); ax[0].legend(fontsize=7)
    ax[1].axhline(0, color="k", lw=0.8); ax[1].set_ylabel("zbary sim - data [layer]")
    ax[2].axhline(1, color="k", lw=0.8); ax[2].set_ylabel("front (L0-2) energy fraction, sim / data")
    for a in ax[:3]:
        a.set_xlabel("extra Al upstream [X0]")
    ax[3].axhline(1, color="k", lw=0.8); ax[3].set_xlabel("layer"); ax[3].set_ylabel("f_l sim / data")
    ax[3].set_title("dashed: v7, solid: Al scan point nearest t_min", fontsize=9); ax[3].legend(fontsize=7)
    ax[3].set_ylim(0.6, 1.3)
    fig.tight_layout(); fig.savefig(f"{args.scan}/upstream_fit.png", dpi=110); plt.close(fig)
    open(f"{args.scan}/upstream_fit.txt", "w").write("\n".join(lines) + "\n")
    print("\n".join(lines))


if __name__ == "__main__":
    main()

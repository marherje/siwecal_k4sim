#!/usr/bin/env python
"""
sim_tuning 1b: where along the beam line is the missing upstream material? Same amount (0.20 X0), different
distances, and an H2-like inventory, against data.

    python analysis/upstream_position.py <outdir> --sim "Al -2.2 m=<dir>" --sim "Al -30 m=<dir>" ... \
        [--control "up0p20=<dir>,up0p20far=<dir>"] [--groups ...]

Per group and sample (digitised valtrees, shower-like events, as upstream_profile.py / tuning_lateral.py):
  R90, f(r>22 mm)   per layer, the lateral width about the event barycentre
  profile chi2      normalised longitudinal profile vs data (sim variance at data statistics)
  hits/ev, E/hit    the digitised totals
--control a=dir,b=dir: the same geometry with the beam born at -2.65 m (a) and far upstream (b); must agree if the
digitiser's TimeOffsetNs undoes the extra flight time.
Writes position.txt, position_width.png (R90 and f>22 per layer, sim - data, one panel per group),
position_summary.png (chi2, width of L0-2 and hits/ev per sample and group) and, with --control, control.txt.
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
from analysis.tuning_lateral import lateral, r90, RB  # noqa: E402

NL = st.NL
NSHOW = 9
COLORS = ["#8c8c8c", "#eda100", "#1baf7a", "#c0392b", "#4a3aa7", "#17becf"]


def summary(s):
    L = lateral(s)
    R = np.array([r90(L["prof"][l]) for l in range(NL)])
    f22 = np.array([L["prof"][l][RB[1:] > 22].sum() / max(L["prof"][l].sum(), 1e-12) for l in range(NL)])
    nev = int(s["nev"])
    # normalised longitudinal profile per event
    ev, sl, e = s["ev"], s["slab"].astype(int), s["e"].astype(float)
    epl = np.zeros((nev, NL))
    np.add.at(epl, (ev, sl), e)
    tot = epl.sum(1)
    ok = tot > 0
    f = epl[ok] / tot[ok, None]
    return {"R90": R, "f22": f22, "f": f.mean(0), "fvar": f.var(0) / ok.sum(), "n": ok.sum(),
            "hits": len(e) / nev, "eph": e.mean(), "E": e.sum() / nev}


def chi2(r, d):
    return float(((r["f"] - d["f"]) ** 2 / (r["fvar"] * r["n"] / d["n"] + d["fvar"])).sum() / NL)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("outdir")
    ap.add_argument("--sim", action="append", required=True)
    ap.add_argument("--control", default=None)
    ap.add_argument("--groups", default="th210_e20_26_-52,th230_E4_e10,th230_E4_e20,th230_E4_e74")
    args = ap.parse_args()
    os.makedirs(args.outdir, exist_ok=True)
    groups = {g["name"]: g for g in st.load_groups()}
    names = args.groups.split(",")
    sims = [s.split("=", 1) for s in args.sim]
    tag = lambda lab: "pos_" + "".join(c if c.isalnum() else "_" for c in lab)

    res = {}
    lines = ["group sample hits/ev E/hit chi2 | R90 L0..L8 | f22 L0..L8"]
    for n in names:
        d = st.load(st.data_paths(groups[n]), cache=f"{st.OUT}/cache/data_{n}.npz")
        res[n] = {"data": summary(d)}
        for lab, path in sims:
            p = st.sim_paths(n, path)
            if p:
                res[n][lab] = summary(st.load(p, cache=f"{st.OUT}/cache/{tag(lab)}_{n}.npz"))
        for lab, r in res[n].items():
            c = 0.0 if lab == "data" else chi2(r, res[n]["data"])
            r["chi2"] = c
            lines.append(f"{n} {lab:<14} {r['hits']:7.1f} {r['eph']:6.3f} {c:7.0f} | "
                         + " ".join(f"{x:5.1f}" for x in r["R90"][:NSHOW]) + " | "
                         + " ".join(f"{x:.3f}" for x in r["f22"][:NSHOW]))
        print("\n".join(lines[-len(res[n]):]), flush=True)
    open(f"{args.outdir}/position.txt", "w").write("\n".join(lines) + "\n")

    # per-layer width, sim - data
    fig, ax = plt.subplots(2, len(names), figsize=(4.4 * len(names), 6.6), sharex=True, squeeze=False)
    L = np.arange(NSHOW)
    for j, n in enumerate(names):
        D = res[n]["data"]
        for i, (lab, _) in enumerate(sims):
            if lab not in res[n]:
                continue
            r = res[n][lab]
            ax[0, j].plot(L, r["R90"][:NSHOW] - D["R90"][:NSHOW], "-o", ms=3, color=COLORS[i % len(COLORS)], label=lab)
            ax[1, j].plot(L, r["f22"][:NSHOW] - D["f22"][:NSHOW], "-o", ms=3, color=COLORS[i % len(COLORS)])
        for a in ax[:, j]:
            a.axhline(0, color="k", lw=0.8)
        ax[0, j].set_title(n, fontsize=10); ax[1, j].set_xlabel("layer")
    ax[0, 0].set_ylabel("R90 sim - data [mm]"); ax[1, 0].set_ylabel("f(r > 22 mm) sim - data")
    ax[0, 0].legend(fontsize=8)
    fig.suptitle("Lateral width of the first layers: 0.20 X0 at different distances, and the H2-like line")
    fig.tight_layout(); fig.savefig(f"{args.outdir}/position_width.png", dpi=110); plt.close(fig)

    # summary per sample
    fig, ax = plt.subplots(1, 3, figsize=(15, 4.2))
    x = np.arange(len(sims))
    for j, n in enumerate(names):
        D = res[n]["data"]
        ok = [i for i, (lab, _) in enumerate(sims) if lab in res[n]]
        col = plt.cm.viridis(j / max(len(names) - 1, 1))
        ax[0].plot(x[ok], [res[n][sims[i][0]]["chi2"] for i in ok], "-o", color=col, label=n)
        ax[1].plot(x[ok], [res[n][sims[i][0]]["R90"][:3].mean() - D["R90"][:3].mean() for i in ok], "-o", color=col)
        ax[2].plot(x[ok], [res[n][sims[i][0]]["hits"] / D["hits"] for i in ok], "-o", color=col)
    for a in ax:
        a.set_xticks(x); a.set_xticklabels([lab for lab, _ in sims], rotation=25, ha="right", fontsize=8)
    ax[0].set_ylabel("longitudinal profile chi2 / ndf"); ax[0].legend(fontsize=8)
    ax[1].axhline(0, color="k", lw=0.8); ax[1].set_ylabel("mean R90 of L0-2, sim - data [mm]")
    ax[2].axhline(1, color="k", lw=0.8); ax[2].set_ylabel("hits / event, sim / data")
    fig.tight_layout(); fig.savefig(f"{args.outdir}/position_summary.png", dpi=110); plt.close(fig)

    if args.control:
        (la, pa), (lb, pb) = [c.split("=", 1) for c in args.control.split(",")]
        cl = [f"control: {la} (beam at -2.65 m) vs {lb} (beam far upstream, TimeOffsetNs)",
              "group  hits/ev a b  ratio | E/hit a b ratio | R90(L0-2) a b"]
        for n in names:
            ra = res[n].get(la) or summary(st.load(st.sim_paths(n, pa), cache=f"{st.OUT}/cache/{tag(la)}_{n}.npz"))
            rb = res[n].get(lb) or summary(st.load(st.sim_paths(n, pb), cache=f"{st.OUT}/cache/{tag(lb)}_{n}.npz"))
            cl.append(f"{n} {ra['hits']:.1f} {rb['hits']:.1f} {rb['hits'] / ra['hits']:.4f} | {ra['eph']:.3f} "
                      f"{rb['eph']:.3f} {rb['eph'] / ra['eph']:.4f} | {ra['R90'][:3].mean():.1f} {rb['R90'][:3].mean():.1f}")
        open(f"{args.outdir}/control.txt", "w").write("\n".join(cl) + "\n")
        print("\n".join(cl))


if __name__ == "__main__":
    main()

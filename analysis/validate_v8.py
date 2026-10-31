#!/usr/bin/env python
"""
Validation of the new default simulation (2026-10-26: H2-like upstream line in SND_compact_beamline.xml, beam at
-65.5 m with RealDigitizer.TimeOffsetNs, AdcDigitizer.OccupancySlope = -1.5e-3) against data and final_v7.

    python analysis/validate_v8.py <outdir> [--new <campaign dir>] [--nmax 40000]

Per beam point of the validation mini-campaign (th210 central scan 7.5-99 GeV, th230 E4 7.5-74 GeV, th230 E5 99 GeV),
shower-like events (nhit_chan > 0.5 x p90, sim_tuning.load):
  hits/ev, E/hit           mean nhit and sum E / sum nhit of the selected hits (valtree), sim / data
  profile chi2             normalised longitudinal profile vs data (sim variance at data statistics)
  R90 L0-2                 lateral width of the first layers about the barycentre, sim - data
  chip-multiplicity chi2   hits-per-chip spectrum vs data (tuning_occupancy.mult_chi2)
Writes validation.txt, validation.npz, validation_summary.png (the four figures of merit against E, v7 vs new),
validation_trend.png (nhits and E/hit against E, data / v7 / new) and validation_profiles.png (longitudinal profile
and R90 per layer at three energies). Style (analysis/plot_style.py): data blue points with error bars, the two
digitised simulations (final_v7 light red, new dark red) as lines in the trends and as hatched histograms in the
per-layer profile.
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
from analysis.tuning_lateral import lateral, r90  # noqa: E402
from analysis.tuning_occupancy import mult_chi2  # noqa: E402
from analysis import plot_style as ps  # noqa: E402
from analysis.plot_style import INK2  # noqa: E402

NL = st.NL
NEW = f"{st.OUT}/validation"
POINTS = ["th210_e7.5_17_-50", "th210_e10_19_-51", "th210_e20_26_-52", "th210_e34_27_-53", "th210_e52_24_-51",
          "th210_e74_22_-49", "th210_e99_22_-48", "th230_E4_e7.5", "th230_E4_e10", "th230_E4_e20", "th230_E4_e34",
          "th230_E4_e52", "th230_E4_e74", "th230_E5_e99"]
STYLE = dict(zip(("final_v7", "new"), ps.reds(2)))      # final_v7 light red, new dark red
LAB = {"data": "data", "final_v7": "sim_v7", "new": "new default (H2-like line + k$\\cdot$Q)"}

plt.rcParams.update(ps.RC)


def summary(s):
    nev = int(s["nev"])
    ev, sl, e = s["ev"], s["slab"].astype(int), s["e"].astype(float)
    epl = np.zeros((nev, NL))
    np.add.at(epl, (ev, sl), e)
    tot = epl.sum(1)
    ok = tot > 0
    f = epl[ok] / tot[ok, None]
    nh = np.bincount(ev, minlength=nev)
    prof = lateral(s)["prof"]
    return {"f": f.mean(0), "fvar": f.var(0) / ok.sum(), "n": ok.sum(), "hits": nh.mean(),
            "hits_err": nh.std() / np.sqrt(nev), "eph": e.sum() / len(e), "eph_err": e.std() / np.sqrt(len(e)),
            "R90": np.array([r90(prof[l]) for l in range(NL)])}


def chi2(r, d):
    return float(((r["f"] - d["f"]) ** 2 / (r["fvar"] * r["n"] / d["n"] + d["fvar"])).sum() / NL)


def compute(args, groups):
    """Figures of merit per point from the sim_tuning caches; writes validation.txt and validation.npz."""
    camp = {"final_v7": st.V7, "new": args.new}
    R = {}
    lines = ["point E | hits/ev data v7 new | E/hit data v7 new | prof chi2 v7 new | R90(L0-2) data v7 new | "
             "mult chi2 v7 new"]
    for name in POINTS:
        g = groups[name]
        d = st.load(st.data_paths(g), nmax=args.nmax, cache=f"{st.OUT}/cache/data_{name}.npz")
        r = {"data": summary(d)}
        for tag, c in camp.items():
            p = st.sim_paths(name, c)
            if not p:
                continue
            ctag = args.cache_tag if tag == "new" else tag     # one cache per campaign, not per role
            s = st.load(p, nmax=args.nmax, cache=f"{st.OUT}/cache/val_{ctag}_{name}.npz")
            r[tag] = summary(s)
            r[tag]["chi2"] = chi2(r[tag], r["data"])
            r[tag]["mult"] = mult_chi2(s, d)
        if "new" not in r:
            print(f"{name}: no new sample yet")
            continue
        R[name] = (float(g["E"]), g["set"], r)
        w = {k: v["R90"][:3].mean() for k, v in r.items()}
        lines.append(f"{name} {g['E']} | {r['data']['hits']:.1f} {r['final_v7']['hits']:.1f} {r['new']['hits']:.1f} | "
                     f"{r['data']['eph']:.3f} {r['final_v7']['eph']:.3f} {r['new']['eph']:.3f} | "
                     f"{r['final_v7']['chi2']:.0f} {r['new']['chi2']:.0f} | "
                     f"{w['data']:.1f} {w['final_v7']:.1f} {w['new']:.1f} | "
                     f"{r['final_v7']['mult']:.0f} {r['new']['mult']:.0f}")
        print(lines[-1], flush=True)
    open(f"{args.outdir}/validation.txt", "w").write("\n".join(lines) + "\n")

    np.savez(f"{args.outdir}/validation.npz", **{
        f"{n}__{tag}__{q}": np.atleast_1d(R[n][2][tag][q]) for n in R for tag in ("data", "final_v7", "new")
        for q in ("hits", "hits_err", "eph", "eph_err", "R90", "f", "fvar") + (("chi2", "mult") if tag != "data" else ())})

    return R


def load_npz(path, groups):
    """R from a validation.npz of an earlier pass (--replot): redraw without reading the caches."""
    z = np.load(path)
    R = {}
    for key in z.files:
        n, tag, q = key.split("__")
        v = z[key]
        R.setdefault(n, (float(groups[n]["E"]), groups[n]["set"], {}))[2].setdefault(tag, {})[q] = \
            v if q in ("R90", "f", "fvar") else float(v[0])
    return R


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("outdir")
    ap.add_argument("--new", default=NEW)
    ap.add_argument("--nmax", type=int, default=40000)
    ap.add_argument("--cache-tag", default="new", help="cache name of the --new campaign (one per campaign)")
    ap.add_argument("--label", default=None, help="legend label of the --new campaign")
    ap.add_argument("--replot", default=None, help="validation.npz of an earlier pass: redraw only, no cache reading")
    args = ap.parse_args()
    os.makedirs(args.outdir, exist_ok=True)
    if args.label:
        LAB["new"] = args.label
    groups = {g["name"]: g for g in st.load_groups()}
    R = load_npz(args.replot, groups) if args.replot else compute(args, groups)
    sets = {"th210": [n for n in R if R[n][1] == "th210"], "th230": [n for n in R if R[n][1] == "th230"]}

    # figures of merit vs E
    fig, ax = plt.subplots(2, 2, figsize=(11, 7.2))
    for sname, mk in (("th210", "o"), ("th230", "s")):
        ns = sorted(sets[sname], key=lambda n: R[n][0])
        E = np.array([R[n][0] for n in ns])
        dd = lambda key: np.array([R[n][2]["data"][key] for n in ns])
        for tag in ("final_v7", "new"):
            q = lambda key: np.array([R[n][2][tag][key] for n in ns])
            sty = dict(STYLE[tag], marker=mk, ls="-" if sname == "th210" else "--")
            lab = f"{LAB[tag]}, {sname}"
            ps.sim_line(ax[0, 0], E, q("hits") / dd("hits"), sty, label=lab)
            ps.sim_line(ax[0, 1], E, q("eph") / dd("eph"), sty, label=lab)
            ps.sim_line(ax[1, 0], E, q("chi2"), sty, label=lab)
            ps.sim_line(ax[1, 1], E, [R[n][2][tag]["R90"][:3].mean() - R[n][2]["data"]["R90"][:3].mean() for n in ns],
                        sty, label=lab)
    for a in ax[0]:
        a.axhline(1, color=INK2, lw=0.8)
    ax[1, 1].axhline(0, color=INK2, lw=0.8)
    ps.pad_band(ax[1, 1])
    band = [h for h, l in zip(*ax[1, 1].get_legend_handles_labels()) if "half a pad" in l]
    ax[1, 1].legend(band, [h.get_label() for h in band], fontsize=7.5, loc="lower center")
    ax[0, 0].set_ylabel("hits / event, sim / data")
    ax[0, 1].set_ylabel("E / hit, sim / data")
    ax[1, 0].set_ylabel("longitudinal profile $\\chi^2$ / ndf")
    ax[1, 0].set_yscale("log")
    ax[1, 1].set_ylabel("R90 of layers 0-2, sim $-$ data [mm]")
    for a in ax[1]:
        a.set_xlabel("beam energy [GeV]")
    ax[0, 0].legend(fontsize=7.5)
    fig.tight_layout()
    fig.savefig(f"{args.outdir}/validation_summary.png", dpi=120)
    plt.close(fig)

    # nhits and E/hit vs E
    fig, ax = plt.subplots(2, 2, figsize=(11, 7.2), sharex=True)
    for j, sname in enumerate(("th210", "th230")):
        ns = sorted(sets[sname], key=lambda n: R[n][0])
        E = np.array([R[n][0] for n in ns])
        for tag in ("final_v7", "new"):
            ps.sim_line(ax[0, j], E, [R[n][2][tag]["hits"] for n in ns], STYLE[tag], label=LAB[tag])
            ps.sim_line(ax[1, j], E, [R[n][2][tag]["eph"] for n in ns], STYLE[tag])
        d = lambda key: np.array([R[n][2]["data"][key] for n in ns])
        ps.data_points(ax[0, j], E, d("hits"), d("hits_err"), label=LAB["data"])
        ps.data_points(ax[1, j], E, d("eph"), d("eph_err"), label=None)
        ax[0, j].set_title(f"{sname} ({'central scan' if sname == 'th210' else 'position E4, E5 at 99 GeV'})",
                           fontsize=10)
        ax[1, j].set_xlabel("beam energy [GeV]")
    ax[0, 0].set_ylabel("hits / event")
    ax[1, 0].set_ylabel("E / hit [MIP]")
    ax[0, 0].legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(f"{args.outdir}/validation_trend.png", dpi=120)
    plt.close(fig)

    # profiles at three energies
    show = [n for n in ("th230_E4_e10", "th230_E4_e20", "th230_E4_e74") if n in R]
    fig, ax = plt.subplots(2, len(show), figsize=(4 * len(show), 6.4), squeeze=False)
    L = np.arange(NL)
    edges = np.arange(NL + 1) - 0.5
    for j, n in enumerate(show):
        r = R[n][2]
        for tag in ("final_v7", "new"):
            ps.sim_hist(ax[0, j], r[tag]["f"], edges, STYLE[tag], label=LAB[tag], fill=(tag == "new"))
            ps.sim_line(ax[1, j], L[:9], r[tag]["R90"][:9] - r["data"]["R90"][:9], STYLE[tag])
        ps.data_points(ax[0, j], L, r["data"]["f"], np.sqrt(r["data"]["fvar"]), edges, label=LAB["data"])
        ax[1, j].axhline(0, color=INK2, lw=1)
        ps.pad_band(ax[1, j], label=(j == 0))
        ax[0, j].set_title(f"{n.replace('_', ' ')} GeV".replace(" e", " "), fontsize=10)
        ax[1, j].axhline(0, color=INK2, lw=0.8)
        ax[1, j].set_xlabel("layer")
    ax[0, 0].set_ylabel("fraction of the event energy")
    ax[1, 0].set_ylabel("R90 sim $-$ data [mm]")
    ax[1, 0].legend(fontsize=8, loc="lower right")
    h, l = ax[0, 0].get_legend_handles_labels()   # one row on top: inside a panel it hides the histograms
    fig.legend(h, l, loc="upper center", ncol=len(l), fontsize=8.5)
    fig.tight_layout(rect=[0, 0, 1, 0.94])
    fig.savefig(f"{args.outdir}/validation_profiles.png", dpi=120)
    plt.close(fig)


if __name__ == "__main__":
    main()

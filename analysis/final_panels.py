#!/usr/bin/env python
"""
Shower characterisation panels for the final TB2026 sample: data (Reconstructed_final, fixed SCA pairing, th210
pedestals + MIP table for every run) against the simulation without and with digitisation (final_v5).

One 3 x 5 panel per group (threshold set, beam energy, beam position); runs of the same group are pooled.

    python analysis/final_panels.py <outdir> [--groups th210_e20_scan,...] [--no-cache]

Writes <outdir>/panel_<group>.png, final_panels_summary.txt, resolution_vs_energy.png, and a per-group cache
(<outdir>/cache/<group>.npz) so re-plotting does not re-read EOS.

Inputs (environment overrides, for tests on other campaigns):
    FINAL_RECO   data campaign        (default TB2026-06/Reconstructed_final)
    FINAL_SIM    simulation folder    (default .../adc_vs_tb/final_v5; trees/ and pid/ inside)
The simulation of group <set>_e<E>_<pos> is trees/ecal_<set>_e<E>_<pos>{,_simple}.root and
pid/<set>_e<E>_<pos>_{digi,nodigi}/*.edm4hep.root.

Selection: shower-like events, nhit > 0.5 x p90 of the sample, applied separately to the ecal tree (energy block)
and to the PID file (shape block), as in th210all_figs.py. Calibrated energy: straight line sum_w [MIP] = a E + b
per sample kind, fitted on the Gaussian-core means of the th210 scan at 7.5-74 GeV (99 GeV is left out: hadron
contamination), applied to every group.
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
sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(_HERE)), "siwecal-tb2026"))
from analysis.compare_event_level import gauss_core  # noqa: E402

B = "/eos/experiment/drdcalo/siw-ecal/TB2026-06"
RECO = os.environ.get("FINAL_RECO", f"{B}/Reconstructed_final")
SIM = os.environ.get("FINAL_SIM", f"{B}/Simulation/Processed/adc_vs_tb/final_v5")


def _run(r):
    return f"TB2026CERN_eudaq_run_{int(r[1:]):06d}" if str(r).startswith("e") else f"TB2026CERN_run_{int(r):06d}"


# (set, energy, position label, pooled runs). The position label is the sim sample's (x_y) beam centre in mm,
# the median shower barycentre of the group in Reconstructed_final (final_v5/README.txt).
GROUPS_FILE = os.path.join(_HERE, "final_v5_groups.txt")


def load_groups(path=GROUPS_FILE):
    """final_v5_groups.txt: '<set> <E> <pos> <runs,comma-separated>' per line ('e' prefix = eudaq run)."""
    out = []
    for line in open(path):
        line = line.split("#")[0].split()
        if len(line) < 4:
            continue
        th, E, pos, runs = line[:4]
        out.append({"set": th, "E": E, "pos": pos, "runs": [_run(r) for r in runs.split(",")],
                    "name": f"{th}_e{E}_{pos}"})
    return out


KINDS = ["data", "nodigi", "digi"]
LABEL = {"data": "Data", "nodigi": "Sim (no digitisation)", "digi": "Sim (digitised)"}
COLOR = {"data": "#29a3dc", "nodigi": "#eda100", "digi": "#c0392b"}
ENERGY_VARS = ["nhit_chan", "sum_hg", "sum_energy", "sum_w"]
SHAPE_VARS = ["e_over_nhit", "moliere", "transverse_rms", "bar_r", "zbary", "shower_start", "shower_max",
              "shower_length", "n_layers_hit", "mip_likeness", "fractal_dimension"]
PANELS = [("nhit_chan", "Number of hits"), ("sum_hg", "ADC sum (high gain) [ADC]"), ("sum_energy", "Energy [MIP]"),
          ("gev", "Calibrated energy [GeV]"), ("e_over_nhit", "Energy per hit [MIP]"),
          ("moliere", "Molière radius (90 %) [mm]"), ("transverse_rms", "Transverse RMS [mm]"),
          ("bar_r", "Barycentre radius [mm]"), ("zbary", "Longitudinal barycentre [layer]"),
          ("shower_start", "Shower start [layer]"), ("shower_max", "Shower maximum [layer]"),
          ("shower_length", "Shower length [layers]"), ("n_layers_hit", "Layers hit"),
          ("mip_likeness", "MIP-likeness"), ("fractal_dimension", "Fractal dimension")]
FIT_PANELS = ("nhit_chan", "sum_hg", "sum_energy", "gev")
INTEGER_VARS = ("shower_start", "shower_max", "shower_length", "n_layers_hit")
CAL_E = ["7.5", "10", "20", "34", "52", "74"]


def showerlike(n):
    return n > 0.5 * np.percentile(n, 90)


def sim_paths(g, kind):
    tree = f"{SIM}/trees/ecal_{g['name']}{'_simple' if kind == 'nodigi' else ''}.root"
    pid = glob.glob(f"{SIM}/pid/{g['name']}_{kind}/*.edm4hep.root")
    return tree, (pid[0] if pid else None)


def load_ecal(paths, data):
    br = ["nhit_chan", "sum_hg", "sum_energy"] + (["sum_w_energy"] if data else ["hit_w_energy"])
    cols = {b: [] for b in ENERGY_VARS}
    for p in paths:
        for a in uproot.iterate(f"{p}:ecal", br, step_size="300 MB", library="np"):
            cols["nhit_chan"].append(a["nhit_chan"]); cols["sum_hg"].append(a["sum_hg"])
            cols["sum_energy"].append(a["sum_energy"])
            cols["sum_w"].append(a["sum_w_energy"] if data else np.array([h.sum() for h in a["hit_w_energy"]]))
    c = {b: np.concatenate(v).astype(float) for b, v in cols.items()}
    k = showerlike(c["nhit_chan"])
    return {b: v[k] for b, v in c.items()}


def load_pid(paths):
    from siwecal_common.edm4hep_pid import PidFileReader
    parts = []
    for p in paths:
        r = PidFileReader(p)
        cols = r.scalar_columns()
        parts.append({v: np.asarray(cols[v], float) if v in cols else np.full(len(cols["nhit"]), np.nan)
                      for v in ["nhit"] + SHAPE_VARS})
    c = {v: np.concatenate([p[v] for p in parts]) for v in parts[0]}
    k = showerlike(c["nhit"])
    return {v: x[k] for v, x in c.items()}


def group_samples(g, outdir, use_cache=True):
    cache = f"{outdir}/cache/{g['name']}.npz"
    keys = ENERGY_VARS + SHAPE_VARS
    if use_cache and os.path.exists(cache):
        z = np.load(cache)
        return {k: {v: z[f"{k}__{v}"] for v in keys} for k in KINDS}
    s = {}
    runs = [r for r in g["runs"] if os.path.exists(f"{RECO}/{r}/ecal_{r}.root")]
    if len(runs) != len(g["runs"]):
        print(f"[{g['name']}] WARNING: missing runs {sorted(set(g['runs']) - set(runs))}")
    s["data"] = {**load_ecal([f"{RECO}/{r}/ecal_{r}.root" for r in runs], True),
                 **load_pid([f"{RECO}/{r}/ecal_{r}.edm4hep.root" for r in runs])}
    for kind in ("nodigi", "digi"):
        tree, pid = sim_paths(g, kind)
        s[kind] = {**load_ecal([tree], False), **load_pid([pid])}
    for kind in KINDS:
        s[kind].pop("nhit", None)
    os.makedirs(os.path.dirname(cache), exist_ok=True)
    np.savez(cache, **{f"{k}__{v}": s[k][v] for k in KINDS for v in keys})
    return s


def core_curve(x, bins, mu, sg):
    """The fitted Gaussian core, scaled to the density histogram over its fit window."""
    d, _ = np.histogram(x, bins=bins, density=True)
    ctr = 0.5 * (bins[1:] + bins[:-1])
    w = (ctr > mu - 1.5 * sg) & (ctr < mu + 2 * sg)
    g = np.exp(-0.5 * ((ctr - mu) / sg) ** 2)
    amp = (d[w] * g[w]).sum() / max((g[w] ** 2).sum(), 1e-30)
    xx = np.linspace(mu - 1.5 * sg, mu + 2 * sg, 120)
    return xx, amp * np.exp(-0.5 * ((xx - mu) / sg) ** 2)


def draw_hist(ax, x, bins, kind, label):
    if kind == "data":
        ax.hist(x, bins=bins, density=True, histtype="stepfilled", color=COLOR[kind], alpha=0.22, lw=0)
        ax.hist(x, bins=bins, density=True, histtype="step", color=COLOR[kind], lw=2.0, label=label)
    elif kind == "nodigi":
        ax.hist(x, bins=bins, density=True, histtype="step", color=COLOR[kind], lw=2.0, ls=":", label=label)
    else:
        ax.hist(x, bins=bins, density=True, histtype="stepfilled", facecolor="none", edgecolor=COLOR[kind],
                hatch="//", lw=0, alpha=0.8)
        ax.hist(x, bins=bins, density=True, histtype="step", color=COLOR[kind], lw=1.6, label=label)


def panel(g, s, outpath):
    fig, axes = plt.subplots(3, 5, figsize=(24, 13.5))
    stats = {}
    for ax, (v, xl) in zip(axes.ravel(), PANELS):
        kinds = [k for k in KINDS if not (v == "sum_hg" and k == "nodigi")]
        xs = {k: s[k][v][np.isfinite(s[k][v])] for k in kinds}
        allx = np.concatenate([xs[k] for k in kinds])
        if not len(allx):
            ax.set_axis_off(); continue
        lo, hi = np.percentile(allx, [0.5, 99.5])
        if v in INTEGER_VARS:
            bins = np.arange(np.floor(lo) - 0.5, np.ceil(hi) + 1.5)
        else:
            bins = np.linspace(lo, hi, 55)
        for k in kinds:
            x = xs[k]
            if v in FIT_PANELS:
                mu, sg = gauss_core(x)
                stats[(v, k)] = (mu, sg / mu, len(x))
                fmt = ".1f" if v == "gev" else ".0f"
                lab = f"{LABEL[k]}: $\\mu$ = {mu:{fmt}}, $\\sigma/\\mu$ = {100 * sg / mu:.1f} %, N = {len(x)}"
                draw_hist(ax, x, bins, k, lab)
                if np.isfinite(mu):
                    xx, yy = core_curve(x, bins, mu, sg)
                    ax.plot(xx, yy, color=COLOR[k], lw=1.3, ls="-" if k != "nodigi" else "--")
            else:
                med = float(np.median(x)) if len(x) else np.nan
                stats[(v, k)] = (med, np.nan, len(x))
                draw_hist(ax, x, bins, k, f"{LABEL[k]} (median {med:.2f})")
        if v == "gev":
            ax.axvline(float(g["E"]), color="0.45", lw=1.2, ls="--", label=f"beam energy {g['E']} GeV")
        ax.set_xlabel(xl, fontsize=11)
        ax.grid(alpha=0.25)
        ax.tick_params(labelsize=9)
        if v in FIT_PANELS:
            title = "Gaussian core fit"
            if v == "sum_hg":
                title += " (undigitised sim has no ADC)"
            ax.legend(fontsize=7.6, title=title, title_fontsize=7.8, loc="upper left")
        else:
            ax.legend(fontsize=7.6, loc="best")
    nruns = ", ".join(("eudaq " if "eudaq" in r else "run ") + str(int(r.split("_")[-1])) for r in g["runs"])
    fig.suptitle(f"{g['E']} GeV e$^\\pm$, {g['set']}, beam at ({g['pos'].replace('_', ', ')}) mm  —  data: {nruns}; "
                 f"shower-like events", fontsize=14)
    fig.tight_layout(rect=[0, 0, 1, 0.965])
    fig.savefig(outpath, dpi=100)
    plt.close(fig)
    return stats


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("outdir")
    p.add_argument("--groups", default=None, help="Comma-separated group names (default: all)")
    p.add_argument("--groups-file", default=GROUPS_FILE)
    p.add_argument("--no-cache", action="store_true")
    args = p.parse_args(argv)
    os.makedirs(args.outdir, exist_ok=True)
    groups = load_groups(args.groups_file)
    S = {g["name"]: group_samples(g, args.outdir, not args.no_cache) for g in groups}

    # MIP -> GeV per sample kind on the th210 scan
    scan = {g["E"]: g["name"] for g in groups if g["set"] == "th210" and g["pos"] == SCAN_POS(groups)}
    cal = {}
    for k in KINDS:
        Es = np.array([float(E) for E in CAL_E if E in scan])
        mu = np.array([gauss_core(S[scan[E]][k]["sum_w"])[0] for E in CAL_E if E in scan])
        cal[k] = np.polyfit(Es, mu, 1)
    for s in S.values():
        for k in KINDS:
            s[k]["gev"] = (s[k]["sum_w"] - cal[k][1]) / cal[k][0]

    only = set(args.groups.split(",")) if args.groups else None
    txt = open(f"{args.outdir}/final_panels_summary.txt", "w")

    def out(line=""):
        print(line); txt.write(line + "\n")

    out("MIP -> GeV: sum_w [MIP] = a * E [GeV] + b, fitted on the Gaussian-core mu of the th210 scan 7.5-74 GeV")
    for k in KINDS:
        out(f"   {LABEL[k]:<24} a = {cal[k][0]:8.2f} MIP/GeV   b = {cal[k][1]:8.1f} MIP")
    out("\nPer group: Gaussian-core mu and sigma/mu for the first four variables, median for the others;"
        " ratio = sim / data")
    RES = []
    for g in groups:
        if only and g["name"] not in only:
            continue
        st = panel(g, S[g["name"]], f"{args.outdir}/panel_{g['name']}.png")
        out(f"\n== {g['name']}  ({len(g['runs'])} run(s); N data/nodigi/digi = "
            f"{st[('nhit_chan', 'data')][2]}/{st[('nhit_chan', 'nodigi')][2]}/{st[('nhit_chan', 'digi')][2]})")
        out(f"   {'variable':<18}{'data':>10}{'nodigi':>10}{'digi':>10}{'res data':>10}{'res nodg':>10}"
            f"{'res digi':>10}{'nodigi/data':>13}{'digi/data':>11}")
        for v, _ in PANELS:
            m = {k: st.get((v, k), (np.nan, np.nan, 0)) for k in KINDS}
            r = lambda k: m[k][0] / m["data"][0] if m["data"][0] else np.nan  # noqa: E731
            out(f"   {v:<18}" + "".join(f"{m[k][0]:10.3g}" for k in KINDS) +
                "".join(f"{100 * m[k][1]:9.2f}%" if np.isfinite(m[k][1]) else f"{'':>10}" for k in KINDS) +
                f"{r('nodigi'):13.3f}{r('digi'):11.3f}")
        RES.append((g, {k: st[("gev", k)] for k in KINDS}))
    txt.close()
    resolution_plot(RES, f"{args.outdir}/resolution_vs_energy.png")


def SCAN_POS(groups):
    """Position of the th210 energy scan: the th210 position with the most energies."""
    from collections import Counter
    c = Counter(g["pos"] for g in groups if g["set"] == "th210")
    return c.most_common(1)[0][0]


def resolution_plot(RES, path):
    if not RES:
        return
    sets = sorted({g["set"] for g, _ in RES})
    fig, axes = plt.subplots(1, 2, figsize=(14, 5.6))
    marker = {"th210": "o", "th220": "s", "th230": "^"}
    for g, st in RES:
        E = float(g["E"])
        for k in KINDS:
            mu, res, n = st[k]
            kw = dict(color=COLOR[k], marker=marker.get(g["set"], "o"), ms=8, ls="none",
                      mfc=COLOR[k] if k != "nodigi" else "white", mew=1.6)
            axes[0].plot(E, 100 * res, **kw)
            axes[1].plot(E, mu / E, **kw)
    # legend: kinds by colour, sets by marker
    from matplotlib.lines import Line2D
    h = [Line2D([], [], color=COLOR[k], marker="o", ls="none", ms=8, mfc=COLOR[k] if k != "nodigi" else "white",
                mew=1.6, label=LABEL[k]) for k in KINDS]
    h += [Line2D([], [], color="0.3", marker=marker.get(s, "o"), ls="none", ms=8, label=s) for s in sets]
    for ax, yl in zip(axes, ["Energy resolution $\\sigma/\\mu$ [%]", "$\\mu$ / beam energy (linearity)"]):
        ax.set_xlabel("Beam energy [GeV]"); ax.set_ylabel(yl); ax.grid(alpha=0.25)
    axes[0].legend(handles=h, fontsize=9)
    axes[1].axhline(1, color="0.5", lw=1, ls="--")
    fig.suptitle("Calibrated energy (th210 scan MIP → GeV line per sample), every group; shower-like events")
    fig.tight_layout(rect=[0, 0, 1, 0.94]); fig.savefig(path, dpi=110); plt.close(fig)


if __name__ == "__main__":
    sys.exit(main())

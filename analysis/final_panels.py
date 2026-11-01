#!/usr/bin/env python
"""
Shower characterisation panels for the final TB2026 sample: data (Reconstruction, fixed SCA pairing, th210
pedestals + MIP table for every run) against the simulation without and with digitisation (final_v7).

One 3 x 5 panel per group (threshold set, beam energy, beam position); runs of the same group are pooled.

    python analysis/final_panels.py <outdir> [--groups th210_e20_26_-52,th230_E4_e20,...] [--no-cache]
        [--border-cut MM --cal-file <nocut_outdir>/mip_gev_line.txt]

Writes <outdir>/panel_<group>.png, final_panels_summary.txt, resolution_vs_energy.png, mip_gev_line.txt (the fitted
MIP -> GeV line), and a per-group cache (<outdir>/cache_<campaign>[_b<MM>]/<group>.npz) so re-plotting does not
re-read EOS.

--border-cut MM keeps, among the shower-like events, those with barxy_min_border_distance > MM (the barycentre's
distance to the nearest active edge, needs valtrees that carry it), and writes border_survival.txt. The p90 of the
shower-like cut is still taken on the whole sample. The MIP -> GeV calibration applies to every event, so it is not
refitted on the cut sample: --border-cut requires --cal-file, the mip_gev_line.txt of a pass without the cut.

Inputs (environment overrides, for tests on other campaigns):
    FINAL_RECO   data campaign        (default TB2026-06/Reconstruction)
    FINAL_SIM    simulation folder    (default .../Simulation/Processed/final_v7; trees/ and pid/ inside)
The simulation of group <set>_e<E>_<pos> is pid/<set>_e<E>_<pos>_{digi,nodigi}/*.valtree.root; the data
<run>/ecal_<run>.valtree.root.

Inputs are the valtrees (gaudi_jobs/pid_to_valtree.py) of data and simulation: the PID stage's 0.5 MIP hit cut applies to
every variable, identically in both. Selection: shower-like events, nhit > 0.5 x p90 of the sample. Calibrated energy: straight line sum_w [MIP] = a E + b
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
RECO = os.environ.get("FINAL_RECO", f"{B}/Reconstruction")
SIM = os.environ.get("FINAL_SIM", f"{B}/Simulation/Processed/final_v7")


def _run(r):
    return f"TB2026CERN_eudaq_run_{int(r[1:]):06d}" if str(r).startswith("e") else f"TB2026CERN_run_{int(r):06d}"


# One group per simulated beam point (analysis/beam_points.txt): name, set, table label, energy, beam centre (the
# median shower barycentre of the group in Reconstruction), beam sigma and pooled runs.
GROUPS_FILE = os.environ.get("FINAL_GROUPS", os.path.join(_HERE, "beam_points.txt"))


def load_groups(path=GROUPS_FILE):
    """beam_points.txt: '<name> <set> <label> <E> <x> <y> <sx> <sy> <runs,comma-separated>' ('e' prefix = eudaq
    run, label '-' = none). Older files: final_v5_groups.txt '<set> <E> <x_y> <runs>', group <set>_e<E>_<x_y>;
    final_v6_groups.txt '<set> <label> <E> <x_y> <runs>', group <set>_<label>_e<E>."""
    out = []
    for line in open(path):
        line = line.split("#")[0].split()
        if len(line) < 4:
            continue
        if len(line) >= 9:
            name, th, label, E, x, y, _, _, runs = line[:9]
            pos, label = f"{x}_{y}", None if label == "-" else label
        elif len(line) >= 5:
            th, label, E, pos, runs = line[:5]
            name = f"{th}_{label}_e{E}"
        else:
            th, E, pos, runs = line[:4]
            label, name = None, f"{th}_e{E}_{pos}"
        out.append({"set": th, "E": E, "pos": pos, "label": label, "runs": [_run(r) for r in runs.split(",")],
                    "name": name})
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


def sim_valtree(g, kind):
    v = glob.glob(f"{SIM}/pid/{g['name']}_{kind}/*.valtree.root")
    return v[0] if v else None


VT_BRANCHES = ["nhit_chan", "sum_hg", "sum_energy", "weighte"] + SHAPE_VARS


BORDER_VAR = "barxy_min_border_distance"


def load_valtrees(paths, border_cut=None):
    """Every variable from the valtrees (pid_to_valtree.py): the PID stage's hit cut (0.5 MIP) applies to all of
    them, identically in data and simulation. sum_w = weighte (tungsten-weighted energy). With border_cut, the
    shower-like events are further required to have barxy_min_border_distance > border_cut [mm]; n_showerlike
    keeps the count before that cut."""
    branches = VT_BRANCHES + ([BORDER_VAR] if border_cut is not None else [])
    cols = {b: [] for b in branches}
    for p in paths:
        for a in uproot.iterate(f"{p}:ecal", branches, step_size="400 MB", library="np"):
            for b in branches:
                cols[b].append(a[b])
    c = {b: np.concatenate(v).astype(float) for b, v in cols.items()}
    c["sum_w"] = c.pop("weighte")
    k = showerlike(c["nhit_chan"])
    n_showerlike = int(k.sum())
    if border_cut is not None:
        k &= c.pop(BORDER_VAR) > border_cut     # NaN (no barycentre) fails
    out = {b: v[k] for b, v in c.items()}
    out["n_showerlike"] = np.array(n_showerlike)
    return out


def cache_dir(outdir, border_cut=None):
    """Per campaign and border cut, so samples of different campaigns or cuts never mix."""
    tag = os.path.basename(os.path.normpath(SIM)) + (f"_b{border_cut:g}" if border_cut is not None else "")
    return f"{outdir}/cache_{tag}"


def group_samples(g, outdir, use_cache=True, border_cut=None):
    cache = f"{cache_dir(outdir, border_cut)}/{g['name']}.npz"
    keys = ENERGY_VARS + SHAPE_VARS + ["n_showerlike"]
    if use_cache and os.path.exists(cache):
        z = np.load(cache)
        return {k: {v: z[f"{k}__{v}"] for v in keys} for k in KINDS}
    runs = [r for r in g["runs"] if os.path.exists(f"{RECO}/{r}/ecal_{r}.valtree.root")]
    if len(runs) != len(g["runs"]):
        print(f"[{g['name']}] WARNING: missing runs {sorted(set(g['runs']) - set(runs))}")
    s = {"data": load_valtrees([f"{RECO}/{r}/ecal_{r}.valtree.root" for r in runs], border_cut)}
    for kind in ("nodigi", "digi"):
        s[kind] = load_valtrees([sim_valtree(g, kind)], border_cut)
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
    """Simulations as histograms (dotted / hatched) with their statistical error bars, data as solid points with
    Poisson error bars and a horizontal bar across the bin; all unit area over the drawn range."""
    counts, edges = np.histogram(x, bins=bins)
    norm = max(counts.sum(), 1) * np.diff(edges)
    ctr = 0.5 * (edges[1:] + edges[:-1])
    nz = counts > 0
    if kind == "data":
        ax.errorbar(ctr[nz], (counts / norm)[nz], yerr=(np.sqrt(counts) / norm)[nz],
                    xerr=(0.5 * np.diff(edges))[nz], fmt="o", color=COLOR[kind],
                    ms=6.5, mew=0, elinewidth=1.4, capsize=0, label=label, zorder=5)
        return
    if kind == "nodigi":
        ax.hist(x, bins=bins, density=True, histtype="step", color=COLOR[kind], lw=2.0, ls=":", label=label)
    else:
        ax.hist(x, bins=bins, density=True, histtype="stepfilled", facecolor="none", edgecolor=COLOR[kind],
                hatch="//", lw=0, alpha=0.8)
        ax.hist(x, bins=bins, density=True, histtype="step", color=COLOR[kind], lw=1.6, label=label)
    # statistical error per bin, the two simulations side by side so the bars do not overlap
    xs = ctr + (-0.15 if kind == "nodigi" else 0.15) * np.diff(edges)
    ax.errorbar(xs[nz], (counts / norm)[nz], yerr=(np.sqrt(counts) / norm)[nz], fmt="none", ecolor=COLOR[kind],
                elinewidth=1.1, capsize=2, zorder=4)


def selection_label(border_cut=None):
    return "shower-like events" + (f", barxy_min_border_distance > {border_cut:g} mm" if border_cut is not None
                                   else "")


# Threshold sets with their own MIP -> GeV line, fitted on (set, scan label); the other sets use the th210 line.
# The response of a threshold set bends differently with energy (th210 data lose ~9 % of MIP/GeV from 7.5 to 74 GeV,
# th230 ~1 %), so a straight line from the th210 scan does not transfer to th230. th220 has no energy scan of its own.
SET_LINES = {"th230": "E2"}


def read_cal_all(path):
    """mip_gev_line.txt: '<kind> <a [MIP/GeV]> <b [MIP]>' (the th210 line, used by default) and, for a threshold set
    with its own line, '<set>:<kind> <a> <b>'. Returns {None: {kind: (a, b)}, <set>: {kind: (a, b)}, ...}."""
    cals = {}
    for line in open(path):
        f = line.split("#")[0].split()
        if len(f) == 3:
            st, k = f[0].split(":") if ":" in f[0] else (None, f[0])
            cals.setdefault(st, {})[k] = np.array([float(f[1]), float(f[2])])
    for st, cal in cals.items():
        missing = set(KINDS) - set(cal)
        if missing:
            raise SystemExit(f"{path}: no MIP -> GeV line for {st or 'default'} {sorted(missing)}")
    if None not in cals:
        raise SystemExit(f"{path}: no default (th210) MIP -> GeV line")
    return cals


def read_cal(path, threshold_set=None):
    """The MIP -> GeV line of one threshold set (its own if the file has one, else the th210 line)."""
    cals = read_cal_all(path)
    return cals.get(threshold_set, cals[None])


def write_cal(path, cals, source):
    with open(path, "w") as f:
        f.write(f"# MIP -> GeV line, sum_w [MIP] = a * E [GeV] + b, Gaussian-core mu, {CAL_E[0]}-{CAL_E[-1]} GeV; "
                f"{source}\n# default: th210 scan; '<set>:<kind>': that threshold set's own line "
                + ", ".join(f"{st} on {lab}" for st, lab in SET_LINES.items()) + "\n# kind a b\n")
        for st in [None] + [s for s in cals if s is not None]:
            for k in KINDS:
                f.write(f"{(st + ':') if st else ''}{k} {cals[st][k][0]:.6g} {cals[st][k][1]:.6g}\n")


def panel(g, s, outpath, border_cut=None):
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
                 f"{selection_label(border_cut)}", fontsize=14)
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
    p.add_argument("--border-cut", type=float, default=None,
                   help="keep barxy_min_border_distance > MM among the shower-like events")
    p.add_argument("--cal-file", default=None,
                   help="MIP -> GeV line (mip_gev_line.txt of a pass without --border-cut) instead of fitting it")
    args = p.parse_args(argv)
    if args.border_cut is not None and not args.cal_file:
        p.error("--border-cut needs --cal-file: the MIP -> GeV line is fitted without the cut")
    os.makedirs(args.outdir, exist_ok=True)
    groups = load_groups(args.groups_file)
    S = {g["name"]: group_samples(g, args.outdir, not args.no_cache, args.border_cut) for g in groups}

    # MIP -> GeV per sample kind: the th210 scan (default), plus each SET_LINES set on its own scan
    scan = th210_scan(groups)
    scans = {None: scan}
    for st, lab in SET_LINES.items():
        sc = set_scan(groups, st, lab)
        if sc:
            scans[st] = sc
    if args.cal_file:
        cals = read_cal_all(args.cal_file)
        cal_source = f"read from {args.cal_file}"
    else:
        cals = {}
        for st, sc in scans.items():
            cals[st] = {}
            for k in KINDS:
                Es = np.array([float(E) for E in CAL_E if E in sc])
                mu = np.array([gauss_core(S[sc[E]][k]["sum_w"])[0] for E in CAL_E if E in sc])
                cals[st][k] = np.polyfit(Es, mu, 1)
        cal_source = f"fitted on {SIM} and {RECO}"
    cal = cals[None]
    write_cal(f"{args.outdir}/mip_gev_line.txt", cals, cal_source)
    if args.border_cut is not None:
        survival_table(groups, S, args.border_cut, f"{args.outdir}/border_survival.txt")
    for g in groups:
        c = cals.get(g["set"], cal)
        for k in KINDS:
            S[g["name"]][k]["gev"] = (S[g["name"]][k]["sum_w"] - c[k][1]) / c[k][0]

    only = set(args.groups.split(",")) if args.groups else None
    txt = open(f"{args.outdir}/final_panels_summary.txt", "w")

    def out(line=""):
        print(line); txt.write(line + "\n")

    out(f"Selection: {selection_label(args.border_cut)}; simulation {SIM}; data {RECO}")
    out("MIP -> GeV: sum_w [MIP] = a * E [GeV] + b, Gaussian-core mu 7.5-74 GeV, per threshold set: "
        + ", ".join(f"{st} on its {lab} scan" for st, lab in SET_LINES.items() if st in cals)
        + ", every other set on the th210 scan; " + cal_source)
    for st in [None] + [s for s in cals if s is not None]:
        for k in KINDS:
            out(f"   {(st or 'th210 (default)'):<16}{LABEL[k]:<24} a = {cals[st][k][0]:8.2f} MIP/GeV"
                f"   b = {cals[st][k][1]:8.1f} MIP")
    out("\nPer group: Gaussian-core mu and sigma/mu for the first four variables, median for the others;"
        " ratio = sim / data")
    RES = []
    for g in groups:
        if only and g["name"] not in only:
            continue
        st = panel(g, S[g["name"]], f"{args.outdir}/panel_{g['name']}.png", args.border_cut)
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
    resolution_plot(RES, f"{args.outdir}/resolution_vs_energy.png", args.border_cut)
    calibration_plot(S, groups, scan, cal, f"{args.outdir}/calibration_th210.png", args.border_cut)
    for st, sc in scans.items():
        if st is not None and st in cals:
            calibration_plot(S, groups, sc, cals[st], f"{args.outdir}/calibration_{st}.png", args.border_cut,
                             threshold_set=st, scan_label=SET_LINES[st])


def survival_table(groups, S, border_cut, path, min_data=1000):
    """Fraction of the shower-like events with barxy_min_border_distance > border_cut, per group and kind."""
    with open(path, "w") as f:
        f.write(f"# barxy_min_border_distance > {border_cut:g} mm among shower-like events; "
                f"N_pass / N_showerlike (fraction); '!' = fewer than {min_data} data events left\n")
        f.write(f"# {'group':<28}" + "".join(f"{k:>24}" for k in KINDS) + "\n")
        for g in groups:
            s = S[g["name"]]
            cells = []
            for k in KINDS:
                n0, n1 = int(s[k]["n_showerlike"]), len(s[k]["nhit_chan"])
                cells.append(f"{n1:>9d}/{n0:<7d}({n1 / max(n0, 1):5.3f})")
            flag = " !" if len(s["data"]["nhit_chan"]) < min_data else ""
            f.write(f"  {g['name']:<28}" + "".join(f"{c:>24}" for c in cells) + flag + "\n")
    print(open(path).read())


def th210_scan(groups):
    """The th210 energy scan used for the MIP -> GeV line: one th210 group per energy of CAL_E. Selected by
    energy, not by position: every group sits at its own data centroid, so the scan's energies do not share
    one position (selecting by a common position left a single 7.5 GeV point and a meaningless line)."""
    scan = {}
    for g in groups:
        if g["set"] == "th210" and g["E"] in CAL_E:
            if g["E"] in scan:
                raise SystemExit(f"two th210 groups at {g['E']} GeV ({scan[g['E']]}, {g['name']}): "
                                 "the calibration scan needs exactly one per energy")
            scan[g["E"]] = g["name"]
    if len(scan) < 3:
        raise SystemExit(f"th210 calibration scan has only {sorted(scan)} GeV; need at least 3 energies")
    return scan


def set_scan(groups, threshold_set, label):
    """Energy scan of one threshold set at one table position (label, e.g. E2): {E: group name} over CAL_E."""
    return {g["E"]: g["name"] for g in groups
            if g["set"] == threshold_set and g.get("label") == label and g["E"] in CAL_E}


def calibration_plot(S, groups, scan, cal, path, border_cut=None, threshold_set="th210", scan_label=None):
    """The MIP -> GeV line: Gaussian-core mu of sum_w against beam energy on the scan (fitted, 7.5-74 GeV), for th210
    with every th210 99 GeV group drawn on top (not fitted), and the residuals to the line."""
    e99 = [g for g in groups if g["set"] == "th210" and g["E"] == "99"] if threshold_set == "th210" else []
    fig, (ax, axr) = plt.subplots(2, 1, figsize=(8.5, 7.5), sharex=True, gridspec_kw={"height_ratios": [2.2, 1]})
    xx = np.linspace(0, 105, 50)
    for i, k in enumerate(KINDS):
        a, b = cal[k]
        Es = np.array([float(E) for E in CAL_E if E in scan])
        mu = np.array([gauss_core(S[scan[E]][k]["sum_w"])[0] for E in CAL_E if E in scan])
        mu99 = np.array([gauss_core(S[g["name"]][k]["sum_w"])[0] for g in e99])
        dx = (i - 1) * 0.8
        filled = COLOR[k] if k != "nodigi" else "white"
        ax.plot(xx, a * xx + b, color=COLOR[k], lw=1.2, alpha=0.7)
        ax.plot(Es + dx, mu, "o", color=COLOR[k], mfc=filled, mew=1.6, ms=7,
                label=f"{LABEL[k]}: a = {a:.1f} MIP/GeV, b = {b:.0f} MIP")
        ax.plot(np.full(len(mu99), 99.0) + dx, mu99, "D", color=COLOR[k], mfc="none", mew=1.4, ms=7)
        axr.plot(Es + dx, 100 * (mu - (a * Es + b)) / (a * Es + b), "o", color=COLOR[k], mfc=filled, mew=1.6, ms=7)
        axr.plot(np.full(len(mu99), 99.0) + dx, 100 * (mu99 - (a * 99 + b)) / (a * 99 + b), "D", color=COLOR[k],
                 mfc="none", mew=1.4, ms=7)
    if e99:
        ax.plot([], [], "D", color="0.3", mfc="none", label=f"99 GeV, {len(e99)} positions (not fitted)")
    ax.set_ylabel("sum_w, Gaussian-core $\\mu$ [MIP]"); ax.grid(alpha=0.25); ax.legend(fontsize=9, loc="upper left")
    ax.set_title(f"MIP → GeV calibration, {threshold_set} {scan_label or ''} scan: sum_w = a E + b, fitted 7.5–74 GeV" +
                 (f"\npoints: barxy_min_border_distance > {border_cut:g} mm; line without the cut"
                  if border_cut is not None else ""), fontsize=12)
    axr.axhline(0, color="0.5", lw=1, ls="--"); axr.grid(alpha=0.25)
    axr.set_xlabel("Beam energy [GeV]"); axr.set_ylabel("residual to the line [%]")
    fig.tight_layout(); fig.savefig(path, dpi=130); plt.close(fig)


def resolution_plot(RES, path, border_cut=None):
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
    fig.suptitle("Calibrated energy (th210 scan MIP → GeV line per sample), every group; "
                 + selection_label(border_cut))
    fig.tight_layout(rect=[0, 0, 1, 0.94]); fig.savefig(path, dpi=110); plt.close(fig)


if __name__ == "__main__":
    sys.exit(main())

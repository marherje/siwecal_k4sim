#!/usr/bin/env python
"""
Longitudinal shower profiles of the final TB2026 sample against depth in radiation lengths: hits, ADC sum,
energy [MIP] and calibrated energy [GeV] per layer, data (Reconstruction) against the simulation without and
with digitisation (final_v7).

    python analysis/final_profiles.py <outdir> [--groups th210_e20_26_-52,th210_e74_22_-49] [--max-data-events N]
        [--cal-file <final_panels outdir>/mip_gev_line.txt] [--border-cut MM]

Same selection as final_panels.py, rebuilt hit by hit from the ecal trees (the valtrees carry no per-layer ADC):
unmasked hits above 0.5 MIP (the PID stage's hit cut), shower-like events (nhit > 0.5 x p90 of the sample).

Abscissa: hit_X0, the radiation length in front of each layer (0.8 X0 for layer 0, then 1.2 X0 per layer up to
layer 8 and 1.6 X0 from layer 9). Hits, ADC and MIP are what the silicon samples at that depth, drawn per layer and
event. The calibrated energy is the sampling-corrected hit_w_energy divided by the slope a of the MIP->GeV line of
each sample kind (final_panels.py, th210 scan), then by the absorber depth in front of the layer: dE/dX0 [GeV/X0],
whose integral over depth is the shower energy (minus the line's offset b, which is not per layer).

--border-cut MM keeps, among the shower-like events, those with barxy_min_border_distance > MM, taken from the
valtree of the same file (final_panels.py's inputs). The valtree keeps the ecal tree's event order, so the two are
joined row by row, after checking that (run, event, spill, bcid) agree on every row (the simulation repeats event
numbers across chunks, so the key alone is not unique). The data's ecal trees are read from FINAL_RECO_TREES
(default FINAL_RECO), for valtrees regenerated elsewhere. The slope a comes from --cal-file (final_panels.py's
mip_gev_line.txt, fitted without the border cut); without it, the final_v7 slopes below.

Writes <outdir>/profiles_<group>.png and profiles_summary.txt, with a cache per group
(<outdir>/cache_<campaign>[_b<MM>]/).
"""
import argparse
import os
import sys

import numpy as np
import uproot
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(_HERE))
from analysis.final_panels import (RECO, SIM, KINDS, LABEL, COLOR, BORDER_VAR, load_groups, sim_valtree,  # noqa: E402
                                   cache_dir, read_cal, selection_label)

RECO_TREES = os.environ.get("FINAL_RECO_TREES", RECO)

N_LAYERS = 15
HIT_CUT = 0.5
# MIP->GeV slope a [MIP/GeV] per sample kind, th210 scan 7.5-74 GeV (final_panels_summary.txt, final_v7)
SLOPE = {"data": 81.9, "digi": 87.2, "nodigi": 88.7}
MARKER = {"data": "o", "nodigi": "^", "digi": "s"}
LINESTYLE = {"data": "-", "nodigi": ":", "digi": "--"}
QUANTS = [("hits", "Hits per layer"), ("adc", "ADC sum per layer [ADC]"),
          ("mip", "Energy per layer [MIP]"), ("gev", "Calibrated energy, dE/dX$_0$ [GeV / X$_0$]")]
RLIM = (0.6, 1.4)
BRANCHES = ["hit_slab", "hit_energy", "hit_hg", "hit_w_energy", "hit_X0", "hit_ismasked"]
ID_BRANCHES = ["run", "event", "spill", "bcid"]


def sim_tree(g, kind):
    return f"{SIM}/trees/ecal_{g['name']}{'_simple' if kind == 'nodigi' else ''}.root"


def border_distance(tree_path, valtree_path, ids):
    """barxy_min_border_distance of the first len(ids['run']) events of tree_path, from its valtree. Both are in the
    same order, but run_pid_batch.py drops the events without hits, so the valtree is a subsequence of the tree: walk
    the two together, matching (run, event, spill, bcid) row by row; events missing from the valtree get NaN (they
    have no hits and fail the shower-like cut anyway). Every valtree row must find its event, or this stops."""
    n = len(ids["run"])
    v = uproot.open(valtree_path)["ecal"].arrays(ID_BRANCHES + [BORDER_VAR], library="np")
    tkey = np.stack([ids[b] for b in ID_BRANCHES], 1)
    vkey = np.stack([v[b] for b in ID_BRANCHES], 1)
    out = np.full(n, np.nan)
    j = 0
    for i in range(n):
        if j < len(vkey) and (tkey[i] == vkey[j]).all():
            out[i] = v[BORDER_VAR][j]
            j += 1
    full = n == uproot.open(tree_path)["ecal"].num_entries      # whole tree read (no --max-data-events cut)
    if (full and j != len(vkey)) or (n and j == 0):
        raise SystemExit(f"{valtree_path} does not follow {tree_path}: {j} of {len(vkey)} valtree rows matched")
    return out


def layer_sums(paths, max_events=None, valtrees=None, border_cut=None):
    """Per-event, per-layer sums of hits, ADC, MIP and W-weighted MIP after the hit cut; plus the X0 per layer.
    With border_cut, valtrees (one per path) give barxy_min_border_distance for the cut."""
    out = {q: [] for q in ("hits", "adc", "mip", "w")}
    border = []
    x0 = np.full(N_LAYERS, np.nan)
    n_read = 0
    for i, p in enumerate(paths):
        ids = {b: [] for b in ID_BRANCHES}
        for a in uproot.iterate(f"{p}:ecal", BRANCHES + ID_BRANCHES, step_size="300 MB", library="ak"):
            import awkward as ak
            keep = (a["hit_ismasked"] == 0) & (a["hit_energy"] > HIT_CUT)
            sl = a["hit_slab"][keep]
            counts = ak.num(sl)
            ev = np.repeat(np.arange(len(sl)), ak.to_numpy(counts))
            layer = ak.to_numpy(ak.flatten(sl)).astype(int)
            idx = ev * N_LAYERS + layer
            n = len(sl)
            for q, w in (("hits", None), ("adc", a["hit_hg"][keep]), ("mip", a["hit_energy"][keep]),
                         ("w", a["hit_w_energy"][keep])):
                ww = None if w is None else ak.to_numpy(ak.flatten(w)).astype(float)
                out[q].append(np.bincount(idx, weights=ww, minlength=n * N_LAYERS).reshape(n, N_LAYERS))
            xl = ak.to_numpy(ak.flatten(a["hit_X0"][keep]))
            for l in range(N_LAYERS):
                if np.isnan(x0[l]) and np.any(layer == l):
                    x0[l] = float(np.median(xl[layer == l]))
            for b in ID_BRANCHES:
                ids[b].append(ak.to_numpy(a[b]))
            n_read += n
            if max_events and n_read >= max_events:
                break
        if border_cut is not None:
            border.append(border_distance(p, valtrees[i], {b: np.concatenate(v) for b, v in ids.items()}))
        if max_events and n_read >= max_events:
            break
    s = {q: np.concatenate(v)[:max_events] if max_events else np.concatenate(v) for q, v in out.items()}
    nhit = s["hits"].sum(axis=1)
    k = nhit > 0.5 * np.percentile(nhit, 90)
    if border_cut is not None:
        k &= np.concatenate(border)[:len(nhit)] > border_cut     # NaN (no barycentre) fails
    return {q: v[k] for q, v in s.items()}, x0


def group_profiles(g, outdir, max_events, use_cache=True, border_cut=None):
    cache = f"{cache_dir(outdir, border_cut)}/profiles_{g['name']}.npz"
    if use_cache and os.path.exists(cache):
        z = np.load(cache)
        return {k: {q: z[f"{k}__{q}"] for q in ("hits", "adc", "mip", "w")} for k in KINDS}, z["x0"]
    runs = [r for r in g["runs"] if os.path.exists(f"{RECO_TREES}/{r}/ecal_{r}.root")]
    s, x0 = {}, None
    s["data"], x0 = layer_sums([f"{RECO_TREES}/{r}/ecal_{r}.root" for r in runs], max_events,
                               [f"{RECO}/{r}/ecal_{r}.valtree.root" for r in runs], border_cut)
    for kind in ("nodigi", "digi"):
        s[kind], x0k = layer_sums([sim_tree(g, kind)], None, [sim_valtree(g, kind)], border_cut)
        if not np.allclose(x0, x0k, equal_nan=True):
            print(f"[{g['name']}] WARNING: X0 per layer differs, data {x0} {kind} {x0k}")
    os.makedirs(os.path.dirname(cache), exist_ok=True)
    np.savez(cache, x0=x0, **{f"{k}__{q}": s[k][q] for k in KINDS for q in s[k]})
    return s, x0


def profile(s, kind, q, dx0, slope=SLOPE):
    """Mean per event and its statistical error, per layer."""
    v = s[kind]["w"] / slope[kind] / dx0 if q == "gev" else s[kind][q]
    return v.mean(axis=0), v.std(axis=0) / np.sqrt(max(len(v), 1))


def draw_stairs(ax, values, edges, kind, label):
    """A simulation as a histogram with one bin per layer: undigitised dotted, digitised hatched."""
    if kind == "nodigi":
        ax.stairs(values, edges, color=COLOR[kind], lw=2.0, ls=":", label=label)
    else:
        ax.stairs(values, edges, fill=True, facecolor="none", edgecolor=COLOR[kind], hatch="//", lw=0, alpha=0.8)
        ax.stairs(values, edges, color=COLOR[kind], lw=1.6, label=label)


def plot_group(g, s, x0, outpath, slope=SLOPE, border_cut=None):
    dx0 = np.diff(np.r_[0.0, x0])
    # one bin per layer on the depth axis, centred on the layer's X0, edges halfway to its neighbours
    mid = 0.5 * (x0[1:] + x0[:-1])
    edges = np.r_[x0[0] - (mid[0] - x0[0]), mid, x0[-1] + (x0[-1] - mid[-1])]
    fig, axes = plt.subplots(2, 4, figsize=(22, 8.6), gridspec_kw={"height_ratios": [3, 1.25]}, sharex=True)
    summary = []
    for col, (q, ylab) in enumerate(QUANTS):
        ax, rax = axes[0, col], axes[1, col]
        kinds = [k for k in KINDS if not (q == "adc" and k == "nodigi")]
        prof = {k: profile(s, k, q, dx0, slope) for k in kinds}
        md, ed = prof["data"]
        for k in kinds:
            m, e = prof[k]
            tot = (m * dx0).sum() if q == "gev" else m.sum()
            unit = {"hits": "hits", "adc": "ADC", "mip": "MIP", "gev": "GeV"}[q]
            fmt = ".1f" if q == "gev" else ".0f"
            lab = f"{LABEL[k]} (total {tot:{fmt}} {unit}, N = {len(s[k]['hits'])})"
            tmax = x0[np.argmax(m)]
            summary.append((g["name"], q, k, tot, tmax))
            if k == "data":   # solid points with their statistical error bars
                ax.errorbar(x0, m, yerr=e, xerr=[x0 - edges[:-1], edges[1:] - x0], fmt="o", color=COLOR[k],
                            ms=8.5, mew=0, elinewidth=1.5, capsize=0, label=lab, zorder=5)
                continue
            # simulations: one bin per layer, as histograms, with their statistical error per layer
            # (the two simulations side by side in the bin so the bars do not overlap)
            draw_stairs(ax, m, edges, k, lab)
            xs = x0 + (-0.12 if k == "nodigi" else 0.12) * (edges[1:] - edges[:-1])
            ax.errorbar(xs, m, yerr=e, fmt="none", ecolor=COLOR[k], elinewidth=1.3, capsize=2.5, zorder=4)
            r = m / md
            rc = np.clip(r, RLIM[0], RLIM[1])
            draw_stairs(rax, rc, edges, k, None)
            re = np.where(m > 0, r * e / np.maximum(m, 1e-12), 0)   # the simulation's own error; the data's is the band
            inside = r <= RLIM[1]
            rax.errorbar(xs[inside], rc[inside], yerr=re[inside], fmt="none", ecolor=COLOR[k], elinewidth=1.2,
                         capsize=2, zorder=4)
            for xi, ri in zip(x0[r > RLIM[1]], r[r > RLIM[1]]):
                rax.annotate(f"{ri:.1f}", (xi, RLIM[1]), xytext=(4, -4), textcoords="offset points",
                             fontsize=8, color="0.3", va="top")
        # the data's statistical uncertainty, as a band around 1
        rel = np.where(md > 0, ed / np.maximum(md, 1e-12), 0)
        rax.stairs(1 + rel, edges, baseline=1 - rel, fill=True, color="0.6", alpha=0.35, lw=0, label="data stat.")
        ax.set_ylabel(ylab, fontsize=11)
        ax.grid(alpha=0.25)
        ax.legend(fontsize=8, loc="upper right")
        ax.set_ylim(bottom=0)
        rax.axhline(1, color="0.5", lw=1, ls="--")
        rax.set_ylim(*RLIM)
        if col == 0:
            rax.legend(fontsize=8, loc="lower right")
        rax.set_ylabel("sim / data", fontsize=10)
        rax.set_xlabel("Depth in front of the layer [X$_0$]", fontsize=11)
        rax.grid(alpha=0.25)
        top = ax.secondary_xaxis("top")
        top.set_xticks(x0[::2])
        top.set_xticklabels([str(l) for l in range(0, N_LAYERS, 2)], fontsize=8.5)
        top.set_xlabel("Layer", fontsize=9.5)
    nruns = ", ".join(("eudaq " if "eudaq" in r else "run ") + str(int(r.split("_")[-1])) for r in g["runs"])
    fig.suptitle(f"Longitudinal profiles, {g['E']} GeV e$^\\pm$, {g['set']}, beam at ({g['pos'].replace('_', ', ')}) mm"
                 f"  —  data: {nruns}; {selection_label(border_cut)}, hits > {HIT_CUT} MIP", fontsize=13.5)
    fig.tight_layout(rect=[0, 0, 1, 0.955])
    fig.savefig(outpath, dpi=100)
    plt.close(fig)
    return summary


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("outdir")
    ap.add_argument("--groups", default="th210_e20_26_-52,th210_e74_22_-49")
    ap.add_argument("--max-data-events", type=int, default=None)
    ap.add_argument("--no-cache", action="store_true")
    ap.add_argument("--cal-file", default=None, help="mip_gev_line.txt of a final_panels.py pass without the cut")
    ap.add_argument("--border-cut", type=float, default=None,
                    help="keep barxy_min_border_distance > MM among the shower-like events")
    a = ap.parse_args(argv)
    if a.border_cut is not None and not a.cal_file:
        ap.error("--border-cut needs --cal-file: the MIP -> GeV line is fitted without the cut")
    slope = {k: ab[0] for k, ab in read_cal(a.cal_file).items()} if a.cal_file else SLOPE
    os.makedirs(a.outdir, exist_ok=True)
    want = a.groups.split(",")
    groups = [g for g in load_groups() if g["name"] in want]
    rows = []
    for g in groups:
        s, x0 = group_profiles(g, a.outdir, a.max_data_events, not a.no_cache, a.border_cut)
        rows += plot_group(g, s, x0, f"{a.outdir}/profiles_{g['name']}.png", slope, a.border_cut)
        print(f"[{g['name']}] done")
    with open(f"{a.outdir}/profiles_summary.txt", "w") as f:
        f.write("# group quantity kind total_per_event depth_of_max[X0]\n")
        for r in rows:
            f.write(f"{r[0]} {r[1]} {r[2]} {r[3]:.4g} {r[4]:.2f}\n")
    print(open(f"{a.outdir}/profiles_summary.txt").read())


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""
Event-level comparison of a digitised simulation with a test-beam run, on the
same ecal tree, plus the resolution of both.

Two figures and one text summary per point:

  event_overview_<tag>.png   sum_hg [ADC], sum_energy [MIP], hits, slabs hit,
                             and under each its per-layer profile
  shower_shape_<tag>.png     per-hit energy [MIP], radial energy profile, shower
                             depth (energy-weighted layer), hottest cell per event,
                             barycentre x / y
  summary_<tag>.txt          mean, peak, Gaussian core (mu, sigma, sigma/mu) of the
                             three event sums, data and simulation, and the ratios

The resolution is the Gaussian core: fit iterated in [mu - 1.5 sigma, mu + 2 sigma]
around the peak, so the data's low tail (beam contamination, split events) does
not enter sigma; the raw sigma/mu of the full distribution is printed next to it.

Usage
-----
    python -m analysis.compare_event_level --data ecal_TB2026CERN_run_000013.root \\
        --sim ecal_th230_e52.root --energy 52 --tag th230_e52 --title "th230 - e- 52 GeV"
"""

from __future__ import annotations

import argparse
import os
import sys

import numpy as np

C_DATA = "#2a78d6"
C_SIM = "#d2551f"
N_LAYERS = 15


def load(path, max_events=None):
    import uproot

    tree = uproot.open(path)["ecal"]
    branches = ["nhit_chan", "nhit_slab", "sum_hg", "sum_energy", "hit_slab",
                "hit_energy", "hit_hg", "hit_x", "hit_y", "hit_ismasked", "hit_w_energy", "hit_X0"]
    arr = tree.arrays(branches, entry_stop=max_events, library="np")
    # The W-weighted energy (energy x W[slab]/X0, the sampling-corrected sum the
    # test-beam builder writes as sum_w_energy): the simulated tree only carries
    # the per-hit branch, so the event sum is rebuilt here on both sides the
    # same way, from unmasked hits.
    arr["sum_w_energy"] = np.array([float(w[m == 0].sum())
                                    for w, m in zip(arr["hit_w_energy"], arr["hit_ismasked"])])
    return arr


def apply_mip_cut(arr, mip_cut):
    """Drop the hits below `mip_cut` [MIP] and rebuild the event sums from the
    unmasked hits that are left. In place; returns arr."""
    hit_keys = [k for k in arr if k.startswith("hit_")]
    for i in range(len(arr["hit_slab"])):
        keep = arr["hit_energy"][i] > mip_cut
        for k in hit_keys:
            arr[k][i] = arr[k][i][keep]
        m = arr["hit_ismasked"][i] == 0
        arr["nhit_chan"][i] = int(m.sum())
        arr["sum_hg"][i] = float(arr["hit_hg"][i][m].sum())
        arr["sum_energy"][i] = float(arr["hit_energy"][i][m].sum())
        arr["sum_w_energy"][i] = float(arr["hit_w_energy"][i][m].sum())
    return arr


def load_undigitised(path, mip_table=None, mip_cut=0.0, max_events=None, adc_per_mip=None):
    """The simple chain: hit_energy [MIP] from GeV with the per-layer MIP values,
    every cell kept (BasicDigitizer Threshold 0), no trigger, no noise, no ADC.
    Cells below `mip_cut` are dropped (the detector never sees them) and the
    event sums are rebuilt from what is left. The ADC-based observables are
    energy x MPV(slab, chip, channel) from the data's own MIP table, the
    linearised convention of compare_adc_data_sim, so the sample can sit on the
    ADC panels next to the data and the digitised simulation."""
    import uproot
    from analysis.compare_adc_data_sim import mpv_of

    tree = uproot.open(path)["ecal"]
    branches = ["nhit_chan", "nhit_slab", "sum_hg", "sum_energy", "hit_slab", "hit_chip", "hit_chan",
                "hit_energy", "hit_hg", "hit_x", "hit_y", "hit_ismasked", "hit_w_energy", "hit_X0"]
    arr = tree.arrays(branches, entry_stop=max_events, library="np")
    hit_keys = [k for k in branches if k.startswith("hit_")]
    arr["sum_w_energy"] = np.zeros(len(arr["hit_slab"]))
    for i in range(len(arr["hit_slab"])):
        keep = arr["hit_energy"][i] > mip_cut
        for k in hit_keys:
            arr[k][i] = arr[k][i][keep]
        if adc_per_mip is not None:
            # the chip's gain as measured on that set's muons: what an ideal chip
            # (no trigger, no noise, no roll-over) would have output for this energy
            arr["hit_hg"][i] = (arr["hit_energy"][i] * adc_per_mip).astype(np.float32)
        elif mip_table is not None:
            mpv = mpv_of(mip_table, arr["hit_slab"][i].astype(int), arr["hit_chip"][i].astype(int),
                         arr["hit_chan"][i].astype(int))
            arr["hit_hg"][i] = np.where(np.isfinite(mpv), arr["hit_energy"][i] * mpv, 0.0).astype(np.float32)
        m = arr["hit_ismasked"][i] == 0
        arr["nhit_chan"][i] = int(m.sum())
        arr["sum_hg"][i] = float(arr["hit_hg"][i][m].sum())
        arr["sum_energy"][i] = float(arr["hit_energy"][i][m].sum())
        arr["sum_w_energy"][i] = float(arr["hit_w_energy"][i][m].sum())
    return arr


def gauss_core(values, n_iter=5):
    """(mu, sigma) of the Gaussian core around the peak, iterated in an
    asymmetric window [mu - 1.5 sigma, mu + 2 sigma]; nan if too few entries."""
    from scipy.optimize import curve_fit

    v = np.asarray(values, dtype=float)
    v = v[np.isfinite(v) & (v > 0)]
    if len(v) < 50:
        return float("nan"), float("nan")
    # Seed on the bulk: the near-empty events some runs carry (eudaq 286 has a
    # spike at ~800 ADC against a 56k peak) must not be mistaken for the peak.
    bulk = v[v > 0.3 * np.median(v)]
    counts, edges = np.histogram(bulk, bins=60)
    centres = 0.5 * (edges[1:] + edges[:-1])
    mu0 = centres[np.argmax(counts)]
    sigma0 = max(0.5 * np.std(bulk), 1e-3 * mu0)
    mu, sigma = mu0, sigma0
    g = lambda x, a, m, s: a * np.exp(-0.5 * ((x - m) / s) ** 2)   # noqa: E731
    for _ in range(n_iter):
        lo, hi = mu - 1.5 * sigma, mu + 2.0 * sigma
        w = v[(v > lo) & (v < hi)]
        if len(w) < 50:
            break
        counts, edges = np.histogram(w, bins=40)
        centres = 0.5 * (edges[1:] + edges[:-1])
        try:
            (a, m, s), _ = curve_fit(
                g, centres, counts, p0=[counts.max(), mu, sigma],
                bounds=([0, mu0 - 2 * sigma0, 0.05 * sigma0],
                        [np.inf, mu0 + 2 * sigma0, 3 * sigma0]), maxfev=4000)
        except (RuntimeError, ValueError):
            break
        if not np.isfinite(m) or not np.isfinite(s) or s <= 0:
            break
        if abs(m - mu) < 1e-3 * mu and abs(s - sigma) < 1e-3 * sigma:
            mu, sigma = float(m), float(abs(s))
            break
        mu, sigma = float(m), float(abs(s))
    return mu, sigma


def summarise(values):
    v = np.asarray(values, dtype=float)
    v = v[np.isfinite(v)]
    bulk = v[v > 0.3 * np.median(v[v > 0])] if (v > 0).any() else v
    counts, edges = np.histogram(bulk, bins=100)
    peak = 0.5 * (edges[1:] + edges[:-1])[np.argmax(counts)] if counts.sum() else float("nan")
    mu, sigma = gauss_core(v)
    return {"mean": float(v.mean()), "raw_res": float(v.std() / v.mean()) if v.mean() else float("nan"),
            "peak": float(peak), "mu": mu, "sigma": sigma,
            "res": sigma / mu if mu else float("nan")}


def per_layer(arr, weight):
    """Per-event mean of `weight` summed per layer -> array over layers."""
    out = np.zeros(N_LAYERS)
    n = len(arr["hit_slab"])
    for slab, w, m in zip(arr["hit_slab"], weight, arr["hit_ismasked"]):
        keep = m == 0
        out += np.bincount(slab[keep].astype(int), weights=w[keep], minlength=N_LAYERS)[:N_LAYERS]
    return out / max(n, 1)


def shower_shapes(arr):
    """Per-event: depth (energy-weighted layer), hottest cell [MIP], barycentre
    x, y [mm]; and the pooled radial energy profile bins."""
    depth, hottest, bx, by = [], [], [], []
    r_all, e_all = [], []
    for slab, e, x, y, m in zip(arr["hit_slab"], arr["hit_energy"], arr["hit_x"],
                                arr["hit_y"], arr["hit_ismasked"]):
        keep = (m == 0) & (e > 0) & np.isfinite(x) & np.isfinite(y)
        if keep.sum() == 0:
            continue
        e, x, y, slab = e[keep].astype(float), x[keep], y[keep], slab[keep]
        tot = e.sum()
        cx, cy = np.average(x, weights=e), np.average(y, weights=e)
        depth.append(np.average(slab, weights=e))
        hottest.append(e.max())
        bx.append(cx)
        by.append(cy)
        r_all.append(np.hypot(x - cx, y - cy))
        e_all.append(e / tot)
    return (np.array(depth), np.array(hottest), np.array(bx), np.array(by),
            np.concatenate(r_all), np.concatenate(e_all))


def _hist(ax, samples, bins, xlabel, log=False):
    """Unit-area histograms of every sample on one axis, with their core fits.
    `samples` = list of (values, colour, label, stats-or-None, dashed)."""
    width = bins[1] - bins[0]
    x = np.linspace(bins[0], bins[-1], 400)
    for values, colour, label, stats, dashed in samples:
        h, _ = np.histogram(values, bins=bins)
        ax.step(bins[:-1], h / max(h.sum(), 1) / width, where="post", color=colour,
                linewidth=1.6, linestyle="--" if dashed else "-", label=f"{label} ({len(values)} ev)")
        if stats is not None and np.isfinite(stats["mu"]) and np.isfinite(stats["sigma"]):
            g = np.exp(-0.5 * ((x - stats["mu"]) / stats["sigma"]) ** 2)
            g *= 1.0 / (stats["sigma"] * np.sqrt(2 * np.pi))
            ax.plot(x, g, color=colour, linestyle=":", linewidth=1.0, alpha=0.8,
                    label=f"core: μ={stats['mu']:.0f}, σ/μ={stats['res']:.3f}")
    ax.set_xlabel(xlabel)
    ax.set_ylabel("events (unit area)")
    if log:
        ax.set_yscale("log")
    ax.legend(fontsize=7)
    ax.grid(alpha=0.3)


def main(argv=None) -> int:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--data", required=True)
    p.add_argument("--sim", required=True)
    p.add_argument("--energy", type=float, default=None, help="Beam energy [GeV]")
    p.add_argument("--tag", default="cmp")
    p.add_argument("--title", default="")
    p.add_argument("--x-axis", choices=("layer", "x0"), default="layer",
                   help="abscissa of the per-layer profiles: layer number, or the cumulative "
                        "radiation length in front of the layer (hit_X0, from the tree)")
    p.add_argument("--sim-undigi", default=None,
                   help="ecal tree of the UNDIGITISED (simple) chain, drawn as a third curve")
    p.add_argument("--mip-file", default=None,
                   help="data MIP table (layer chip channel mpv); gives the undigitised sample "
                        "its ADC as energy x MPV. Without it its ADC panels stay empty")
    p.add_argument("--undigi-adc-per-mip", type=float, default=None,
                   help="ADC of the undigitised sample as energy x this gain (the set's muon gain, "
                        "ADC/MIP) instead of the data's MIP table")
    p.add_argument("--mip-cut", type=float, default=0.0,
                   help="drop the hits below this energy [MIP] on EVERY sample (data, digitised, "
                        "undigitised) and rebuild the event sums from what is left")
    p.add_argument("--undigi-mip-cut", type=float, default=None,
                   help="cut for the undigitised chain alone (default: --mip-cut, or 0.5 if that is 0: "
                        "its Threshold-0 cells are never seen by the detector)")
    p.add_argument("--outdir", default="plots")
    p.add_argument("--max-data-events", type=int, default=20000)
    args = p.parse_args(argv)

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    os.makedirs(args.outdir, exist_ok=True)
    d = load(args.data, args.max_data_events)
    s = load(args.sim)
    if args.mip_cut > 0:
        d = apply_mip_cut(d, args.mip_cut); s = apply_mip_cut(s, args.mip_cut)
    if args.undigi_mip_cut is None:
        args.undigi_mip_cut = args.mip_cut if args.mip_cut > 0 else 0.5
    u = None
    if args.sim_undigi:
        table = None
        if args.mip_file:
            from analysis.compare_adc_data_sim import load_mip_table
            table = load_mip_table(args.mip_file)
        u = load_undigitised(args.sim_undigi, table, args.undigi_mip_cut, adc_per_mip=args.undigi_adc_per_mip)
    C_UNDIGI = "#17795a"
    U_LABEL = "sim (no digitisation)"

    stats = {}
    for key, label in (("sum_hg", "ADC"), ("sum_energy", "MIP"), ("nhit_chan", "hits"),
                       ("sum_w_energy", "W-weighted")):
        stats[key] = (summarise(d[key]), summarise(s[key])) + ((summarise(u[key]),) if u is not None else ())

    def samples(key):
        out = [(d[key], C_DATA, "data", stats[key][0], False), (s[key], C_SIM, "sim (digi)", stats[key][1], False)]
        if u is not None and not (key == "sum_hg" and args.mip_file is None and args.undigi_adc_per_mip is None):
            out.append((u[key], C_UNDIGI, U_LABEL, stats[key][2], True))
        return out

    # ------------------------------------------------------------ figure 1
    # One column per observable, in the order the reconstruction produces them:
    # ADC -> hits -> energy [MIP] -> W-weighted energy (energy x W/X0 per layer).
    fig, axes = plt.subplots(2, 4, figsize=(19, 8.4))
    lim = lambda k, q=99.5: max([np.percentile(a[k], q) for a in ((d, s) + ((u,) if u is not None else ()))]) * 1.05  # noqa: E731
    _hist(axes[0, 0], samples("sum_hg"), np.linspace(0, lim("sum_hg"), 80),
          "event ADC sum, high gain, pedestal-subtracted")
    _hist(axes[0, 1], samples("nhit_chan"), np.linspace(0, lim("nhit_chan"), 80), "hits per event")
    _hist(axes[0, 2], samples("sum_energy"), np.linspace(0, lim("sum_energy"), 80), "event energy sum [MIP]")
    _hist(axes[0, 3], samples("sum_w_energy"), np.linspace(0, lim("sum_w_energy"), 80),
          "event W-weighted energy sum [MIP x W/X0]")
    # per-layer profiles, one under each event sum, data and sim overlaid; the
    # sim/data total of each is printed in the panel (the per-layer ratio itself
    # is the real_beamline plot's job)
    layers = np.arange(N_LAYERS)
    if args.x_axis == "x0":
        # cumulative X0 of each layer, as the tree carries it per hit (the same on both sides)
        sl_all = np.concatenate(d["hit_slab"][:500]); x0_all = np.concatenate(d["hit_X0"][:500])
        xs = np.array([np.median(x0_all[sl_all == l]) if np.any(sl_all == l) else np.nan for l in layers])
        xlabel = "depth [X0]"
    else:
        xs, xlabel = layers, "layer"
    for ax, weight_key, ylabel in ((axes[1, 0], "hit_hg", "ADC per layer per event"),
                                   (axes[1, 1], None, "hits per layer per event"),
                                   (axes[1, 2], "hit_energy", "energy per layer per event [MIP]"),
                                   (axes[1, 3], "hit_w_energy", "W-weighted energy per layer per event")):
        wd = [np.ones(len(x)) for x in d["hit_slab"]] if weight_key is None else d[weight_key]
        ws = [np.ones(len(x)) for x in s["hit_slab"]] if weight_key is None else s[weight_key]
        pd_, ps_ = per_layer(d, wd), per_layer(s, ws)
        pu_ = None
        if u is not None and not (weight_key == "hit_hg" and args.mip_file is None and args.undigi_adc_per_mip is None):
            wu = [np.ones(len(x)) for x in u["hit_slab"]] if weight_key is None else u[weight_key]
            pu_ = per_layer(u, wu)
        if args.x_axis == "x0":
            # In X0 the sampling cells are not equal (1.2 X0 up to layer 8, 1.6 from
            # layer 9), so a per-layer value drawn as points steps up where the W
            # thickens. Every profile is drawn as a variable-width histogram:
            # height = value / dX0 (a density per X0), AREA of each bin = the
            # per-layer value, and the totals below are the sums of the areas.
            dx0 = np.diff(np.r_[0.0, xs]); edges = np.r_[0.0, xs]
            for prof, colour, name, ls in ((pd_, C_DATA, "data", "-"), (ps_, C_SIM, "sim (digi)", "-"),
                                           (pu_, C_UNDIGI, "sim (no digi)", "--")):
                if prof is not None:
                    ax.stairs(prof / dx0, edges, color=colour, linewidth=1.8, linestyle=ls, label=name)
            ax.set_xlabel(xlabel)
            ax.set_ylabel(ylabel.replace("per layer", "per X0"))
            ax.text(0.97, 0.87, "bin area = per-layer value", transform=ax.transAxes,
                    ha="right", va="top", fontsize=7.5, color="#4a5461")
        else:
            ax.plot(xs, pd_, "o-", color=C_DATA, label="data")
            ax.plot(xs, ps_, "s-", color=C_SIM, label="sim (digi)")
            if pu_ is not None:
                ax.plot(xs, pu_, "^--", color=C_UNDIGI, label="sim (no digi)")
            ax.set_xlabel(xlabel)
            ax.set_ylabel(ylabel)
        ax.legend(fontsize=8)
        ax.grid(alpha=0.3)
        ratio = ps_.sum() / max(pd_.sum(), 1e-9)
        txt = f"sim/data total = {ratio:.3f}"
        if pu_ is not None:
            txt += f"\nno digi / data = {pu_.sum() / max(pd_.sum(), 1e-9):.3f}"
        ax.text(0.97, 0.95, txt, transform=ax.transAxes, ha="right", va="top", fontsize=8.5)
    fig.suptitle(f"Event level — {args.title or args.tag}" + (f" — hits > {args.mip_cut:g} MIP" if args.mip_cut > 0 else ""), fontsize=13)
    fig.tight_layout()
    out1 = os.path.join(args.outdir, f"event_overview_{args.tag}.png")
    fig.savefig(out1, dpi=120)
    plt.close(fig)

    # ------------------------------------------------------------ figure 2
    sd = shower_shapes(d)
    ss = shower_shapes(s)
    fig, axes = plt.subplots(2, 3, figsize=(15, 8.4))
    # per-hit energy spectrum
    ed = np.concatenate([e[m == 0] for e, m in zip(d["hit_energy"], d["hit_ismasked"])])
    es = np.concatenate([e[m == 0] for e, m in zip(s["hit_energy"], s["hit_ismasked"])])
    bins = np.logspace(np.log10(0.3), np.log10(max(ed.max(), es.max(), 1.0)), 60)
    hd, _ = np.histogram(ed, bins=bins)
    hs, _ = np.histogram(es, bins=bins)
    ax = axes[0, 0]
    ax.step(bins[:-1], hd / len(d["nhit_chan"]), where="post", color=C_DATA, label="data")
    ax.step(bins[:-1], hs / len(s["nhit_chan"]), where="post", color=C_SIM, label="sim (digi)")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("hit energy [MIP]")
    ax.set_ylabel("hits per event / bin")
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3)
    # radial energy profile (fraction of the event energy per ring)
    ax = axes[0, 1]
    rbins = np.arange(0, 85, 5)
    for (r, e), colour, label, n_ev in ((sd[4:6], C_DATA, "data", len(sd[0])),
                                        (ss[4:6], C_SIM, "sim (digi)", len(ss[0]))):
        prof, _ = np.histogram(r, bins=rbins, weights=e)
        ax.step(rbins[:-1], prof / max(n_ev, 1), where="post", color=colour, label=label)
    ax.set_yscale("log")
    ax.set_xlabel("distance to the event barycentre [mm]")
    ax.set_ylabel("fraction of the event energy per ring")
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3)
    # depth
    su = shower_shapes(u) if u is not None else None
    two = lambda i: [(sd[i], C_DATA, "data", None, False), (ss[i], C_SIM, "sim (digi)", None, False)] + \
        ([(su[i], C_UNDIGI, "sim (no digi)", None, True)] if su is not None else [])  # noqa: E731
    _hist(axes[0, 2], two(0), np.linspace(0, N_LAYERS, 60), "shower depth: energy-weighted mean layer")
    # hottest cell
    hmax = max(np.percentile(sd[1], 99.5), np.percentile(ss[1], 99.5)) * 1.05
    _hist(axes[1, 0], two(1), np.linspace(0, hmax, 60), "hottest cell per event [MIP]")
    # barycentres
    for ax, i, name in ((axes[1, 1], 2, "x"), (axes[1, 2], 3, "y")):
        lo = min(np.percentile(sd[i], 1), np.percentile(ss[i], 1)) - 10
        hi = max(np.percentile(sd[i], 99), np.percentile(ss[i], 99)) + 10
        _hist(ax, two(i), np.linspace(lo, hi, 60), f"event barycentre {name} [mm]")
    fig.suptitle(f"Shower shape — {args.title or args.tag}", fontsize=13)
    fig.tight_layout()
    out2 = os.path.join(args.outdir, f"shower_shape_{args.tag}.png")
    fig.savefig(out2, dpi=120)
    plt.close(fig)

    # ------------------------------------------------------------ summary
    lines = [f"# {args.title or args.tag}",
             f"#   data: {args.data} ({len(d['nhit_chan'])} events)",
             f"#   sim : {args.sim} ({len(s['nhit_chan'])} events)",
             "" if u is None else f"#   undigitised: {args.sim_undigi} ({len(u['nhit_chan'])} events, cells > {args.undigi_mip_cut:g} MIP)",
             "" if args.energy is None else f"#   beam energy: {args.energy:g} GeV",
             "",
             f"{'quantity':<12}{'side':<6}{'mean':>10}{'peak':>10}{'core mu':>10}"
             f"{'core sigma':>12}{'sigma/mu':>10}{'raw s/m':>9}"]
    for key, label in (("sum_hg", "ADC"), ("nhit_chan", "hits"), ("sum_energy", "MIP"),
                       ("sum_w_energy", "W-weighted")):
        for side, st in (("data", stats[key][0]), ("sim", stats[key][1])) + ((("nodigi", stats[key][2]),) if u is not None else ()):
            lines.append(f"{label:<12}{side:<6}{st['mean']:>10.1f}{st['peak']:>10.1f}"
                         f"{st['mu']:>10.1f}{st['sigma']:>12.1f}{st['res']:>10.4f}"
                         f"{st['raw_res']:>9.3f}")
        sd_, ss_ = stats[key][:2]
        lines.append(f"{'':<12}{'sim/d':<6}{ss_['mean'] / sd_['mean']:>10.3f}"
                     f"{ss_['peak'] / sd_['peak']:>10.3f}{ss_['mu'] / sd_['mu']:>10.3f}"
                     f"{'':>12}{ss_['res'] / sd_['res']:>10.3f}")
        if u is not None:
            su_ = stats[key][2]
            lines.append(f"{'':<12}{'nodg/d':<6}{su_['mean'] / sd_['mean']:>10.3f}"
                         f"{su_['peak'] / sd_['peak']:>10.3f}{su_['mu'] / sd_['mu']:>10.3f}"
                         f"{'':>12}{su_['res'] / sd_['res']:>10.3f}")
    if args.energy:
        lines.append("")
        for key, unit in (("sum_energy", "MIP"), ("sum_hg", "ADC")):
            sd_, ss_ = stats[key][:2]
            lines.append(f"response {unit}/GeV (core mu): data {sd_['mu'] / args.energy:.2f}"
                         f"  sim {ss_['mu'] / args.energy:.2f}")
    lines += ["", f"shower depth [layer]: data {sd[0].mean():.2f}  sim {ss[0].mean():.2f}",
              f"hottest cell [MIP]  : data {np.median(sd[1]):.1f}  sim {np.median(ss[1]):.1f} (median)",
              f"barycentre [mm]     : data ({sd[2].mean():.1f}, {sd[3].mean():.1f})  "
              f"sim ({ss[2].mean():.1f}, {ss[3].mean():.1f})"]
    text = "\n".join(lines)
    print(text)
    with open(os.path.join(args.outdir, f"summary_event_{args.tag}.txt"), "w") as fh:
        fh.write(text + "\n")
    print(f"Written: {out1}\n         {out2}")
    return 0


if __name__ == "__main__":
    sys.exit(main())

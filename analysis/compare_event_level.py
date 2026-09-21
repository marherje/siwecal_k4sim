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
                "hit_energy", "hit_hg", "hit_x", "hit_y", "hit_ismasked", "hit_w_energy"]
    arr = tree.arrays(branches, entry_stop=max_events, library="np")
    # The W-weighted energy (energy x W[slab]/X0, the sampling-corrected sum the
    # test-beam builder writes as sum_w_energy): the simulated tree only carries
    # the per-hit branch, so the event sum is rebuilt here on both sides the
    # same way, from unmasked hits.
    arr["sum_w_energy"] = np.array([float(w[m == 0].sum())
                                    for w, m in zip(arr["hit_w_energy"], arr["hit_ismasked"])])
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


def _hist(ax, d, s, bins, xlabel, log=False, fit=None):
    """Unit-area histograms of data and sim on one axis, with the core fits."""
    hd, _ = np.histogram(d, bins=bins)
    hs, _ = np.histogram(s, bins=bins)
    width = bins[1] - bins[0]
    ax.step(bins[:-1], hd / max(hd.sum(), 1) / width, where="post", color=C_DATA,
            linewidth=1.6, label=f"data ({len(d)} ev)")
    ax.step(bins[:-1], hs / max(hs.sum(), 1) / width, where="post", color=C_SIM,
            linewidth=1.6, label=f"sim (digi) ({len(s)} ev)")
    if fit is not None:
        x = np.linspace(bins[0], bins[-1], 400)
        for stats, colour in ((fit[0], C_DATA), (fit[1], C_SIM)):
            if np.isfinite(stats["mu"]) and np.isfinite(stats["sigma"]):
                g = np.exp(-0.5 * ((x - stats["mu"]) / stats["sigma"]) ** 2)
                g *= 1.0 / (stats["sigma"] * np.sqrt(2 * np.pi))
                # scale to the fraction of entries inside the core
                ax.plot(x, g, color=colour, linestyle="--", linewidth=1.0, alpha=0.8,
                        label=f"core: μ={stats['mu']:.0f}, σ/μ={stats['res']:.3f}")
    ax.set_xlabel(xlabel)
    ax.set_ylabel("events (unit area)")
    if log:
        ax.set_yscale("log")
    ax.legend(fontsize=7.5)
    ax.grid(alpha=0.3)


def main(argv=None) -> int:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--data", required=True)
    p.add_argument("--sim", required=True)
    p.add_argument("--energy", type=float, default=None, help="Beam energy [GeV]")
    p.add_argument("--tag", default="cmp")
    p.add_argument("--title", default="")
    p.add_argument("--outdir", default="plots")
    p.add_argument("--max-data-events", type=int, default=20000)
    args = p.parse_args(argv)

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    os.makedirs(args.outdir, exist_ok=True)
    d = load(args.data, args.max_data_events)
    s = load(args.sim)

    stats = {}
    for key, label in (("sum_hg", "ADC"), ("sum_energy", "MIP"), ("nhit_chan", "hits"),
                       ("sum_w_energy", "W-weighted")):
        stats[key] = (summarise(d[key]), summarise(s[key]))

    # ------------------------------------------------------------ figure 1
    # One column per observable, in the order the reconstruction produces them:
    # ADC -> hits -> energy [MIP] -> W-weighted energy (energy x W/X0 per layer).
    fig, axes = plt.subplots(2, 4, figsize=(19, 8.4))
    lim = lambda k, q=99.5: max(np.percentile(d[k], q), np.percentile(s[k], q)) * 1.05  # noqa: E731
    _hist(axes[0, 0], d["sum_hg"], s["sum_hg"], np.linspace(0, lim("sum_hg"), 80),
          "event ADC sum, high gain, pedestal-subtracted", fit=stats["sum_hg"])
    _hist(axes[0, 1], d["nhit_chan"], s["nhit_chan"], np.linspace(0, lim("nhit_chan"), 80),
          "hits per event", fit=stats["nhit_chan"])
    _hist(axes[0, 2], d["sum_energy"], s["sum_energy"], np.linspace(0, lim("sum_energy"), 80),
          "event energy sum [MIP]", fit=stats["sum_energy"])
    _hist(axes[0, 3], d["sum_w_energy"], s["sum_w_energy"], np.linspace(0, lim("sum_w_energy"), 80),
          "event W-weighted energy sum [MIP x W/X0]", fit=stats["sum_w_energy"])
    # per-layer profiles, one under each event sum, data and sim overlaid; the
    # sim/data total of each is printed in the panel (the per-layer ratio itself
    # is the real_beamline plot's job)
    layers = np.arange(N_LAYERS)
    for ax, weight_key, ylabel in ((axes[1, 0], "hit_hg", "ADC per layer per event"),
                                   (axes[1, 1], None, "hits per layer per event"),
                                   (axes[1, 2], "hit_energy", "energy per layer per event [MIP]"),
                                   (axes[1, 3], "hit_w_energy", "W-weighted energy per layer per event")):
        wd = [np.ones(len(x)) for x in d["hit_slab"]] if weight_key is None else d[weight_key]
        ws = [np.ones(len(x)) for x in s["hit_slab"]] if weight_key is None else s[weight_key]
        pd_, ps_ = per_layer(d, wd), per_layer(s, ws)
        ax.plot(layers, pd_, "o-", color=C_DATA, label="data")
        ax.plot(layers, ps_, "s-", color=C_SIM, label="sim (digi)")
        ax.set_xlabel("layer")
        ax.set_ylabel(ylabel)
        ax.legend(fontsize=8)
        ax.grid(alpha=0.3)
        ratio = ps_.sum() / max(pd_.sum(), 1e-9)
        ax.text(0.97, 0.95, f"sim/data total = {ratio:.3f}", transform=ax.transAxes,
                ha="right", va="top", fontsize=8.5)
    fig.suptitle(f"Event level — {args.title or args.tag}", fontsize=13)
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
    _hist(axes[0, 2], sd[0], ss[0], np.linspace(0, N_LAYERS, 60),
          "shower depth: energy-weighted mean layer")
    # hottest cell
    hmax = max(np.percentile(sd[1], 99.5), np.percentile(ss[1], 99.5)) * 1.05
    _hist(axes[1, 0], sd[1], ss[1], np.linspace(0, hmax, 60), "hottest cell per event [MIP]")
    # barycentres
    for ax, i, name in ((axes[1, 1], 2, "x"), (axes[1, 2], 3, "y")):
        lo = min(np.percentile(sd[i], 1), np.percentile(ss[i], 1)) - 10
        hi = max(np.percentile(sd[i], 99), np.percentile(ss[i], 99)) + 10
        _hist(ax, sd[i], ss[i], np.linspace(lo, hi, 60), f"event barycentre {name} [mm]")
    fig.suptitle(f"Shower shape — {args.title or args.tag}", fontsize=13)
    fig.tight_layout()
    out2 = os.path.join(args.outdir, f"shower_shape_{args.tag}.png")
    fig.savefig(out2, dpi=120)
    plt.close(fig)

    # ------------------------------------------------------------ summary
    lines = [f"# {args.title or args.tag}",
             f"#   data: {args.data} ({len(d['nhit_chan'])} events)",
             f"#   sim : {args.sim} ({len(s['nhit_chan'])} events)",
             "" if args.energy is None else f"#   beam energy: {args.energy:g} GeV",
             "",
             f"{'quantity':<12}{'side':<6}{'mean':>10}{'peak':>10}{'core mu':>10}"
             f"{'core sigma':>12}{'sigma/mu':>10}{'raw s/m':>9}"]
    for key, label in (("sum_hg", "ADC"), ("nhit_chan", "hits"), ("sum_energy", "MIP"),
                       ("sum_w_energy", "W-weighted")):
        for side, st in (("data", stats[key][0]), ("sim", stats[key][1])):
            lines.append(f"{label:<12}{side:<6}{st['mean']:>10.1f}{st['peak']:>10.1f}"
                         f"{st['mu']:>10.1f}{st['sigma']:>12.1f}{st['res']:>10.4f}"
                         f"{st['raw_res']:>9.3f}")
        sd_, ss_ = stats[key]
        lines.append(f"{'':<12}{'sim/d':<6}{ss_['mean'] / sd_['mean']:>10.3f}"
                     f"{ss_['peak'] / sd_['peak']:>10.3f}{ss_['mu'] / sd_['mu']:>10.3f}"
                     f"{'':>12}{ss_['res'] / sd_['res']:>10.3f}")
    if args.energy:
        lines.append("")
        for key, unit in (("sum_energy", "MIP"), ("sum_hg", "ADC")):
            sd_, ss_ = stats[key]
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

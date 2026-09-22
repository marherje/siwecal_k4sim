#!/usr/bin/env python3
"""
Low gain against high gain, simulation against test beam — the dynamic range.

The high-gain preamp saturates; the low gain does not, and the reconstruction
uses it to recover the hits that flattened.  Two things have to be right for that
to work, and both are visible in the same pair of figures:

  the anchor line   adc_low = k * adc_high + c, which EventBuilder::buildHit
                    inverts above SaturationAdc.  In the simulation it is exactly
                    a line, by construction; in the data it is a band with a knee.
  the knee itself   where the high gain stops following the charge.  That is
                    `adc_high_max` in AdcDigitizer, and this script MEASURES it
                    from the data's own hit_hg edge instead of assuming it.

Both sides carry hit_hg and hit_lg in the same branches: the data from
EcalEventBuilder, the simulation from AdcDigitizer's parallel collections through
sim_to_ecal_tree.  So they are compared directly, with no conversion.

Usage
-----
    python -m analysis.compare_gains \\
        --data /eos/.../ecal_TB2026CERN_run_000013.root \\
        --sim  /path/ecal_sim_e52_th230_adc.root \\
        --threshold th230 --tag e52_th230 --outdir plots/
"""

from __future__ import annotations

import argparse
import os
import sys

import numpy as np

# The LG<->HG anchor measured per threshold set, from
# siwecal-tb2026/calibration/MuonCalib_gaudi/anchor/<th>/gain_anchor_<th>.txt.
ANCHOR = {
    "th210": {"k": 0.0925, "c": 2.24, "source_run": "eudaq 166"},
    "th220": {"k": 0.0963, "c": 1.29, "source_run": "run 72"},
    "th230": {"k": 0.0961, "c": 1.60, "source_run": "run 12"},
}

C_DATA = "#29a3dc"
C_SIM = "#c0392b"
C_GRID = "#d8d8d4"
C_INK = "#0b0b0b"
C_INK_SOFT = "#52514e"


def read_hits(path: str, max_events: int | None = None):
    """(hit_hg, hit_lg) of the unmasked hits, both pedestal-subtracted [ADC]."""
    import uproot

    tree = uproot.open(path)["ecal"]
    stop = min(tree.num_entries, max_events) if max_events else tree.num_entries
    arrays = tree.arrays(["hit_hg", "hit_lg", "hit_ismasked"], entry_stop=stop,
                         library="np")
    high = np.concatenate(arrays["hit_hg"]).astype(float)
    low = np.concatenate(arrays["hit_lg"]).astype(float)
    masked = np.concatenate(arrays["hit_ismasked"]).astype(int)
    keep = masked == 0
    if not np.any(low[keep] != 0):
        print(f"[warn] {os.path.basename(path)}: hit_lg is all zeros — a "
              f"simulation made without the ADC model has no low gain")
    return high[keep], low[keep], int(stop)


def roll_over(q, a, n):
    """S(q) = q / (1 + (q/a)^n)^(1/n): the high gain's response to a charge q."""
    return q / np.power(1.0 + np.power(q / a, n), 1.0 / n)


def fit_dynamic_range(high: np.ndarray, low: np.ndarray, k: float, c: float,
                      q_min: float = 300.0) -> dict:
    """Fit the high gain's roll-over from the data's own LG-vs-HG relation.

    The low gain is linear over the whole range, so `(hit_lg - c) / k` is the
    charge in high-gain-equivalent ADC, and `hit_hg` against it IS the response
    curve S(q).  Fitting S(q) = q / (1 + (q/A)^n)^(1/n) to the profile gives the
    two numbers AdcDigitizer needs: the ceiling A and the sharpness n.  This is
    what "anchoring the dynamic range on electrons" means -- only electrons
    reach the knee.

    The profile is the median hit_hg in bins of q (robust to the stray channels
    that sit on other gain ratios), and the fit is weighted by the number of
    entries so the linear region does not dominate the knee.
    """
    from scipy.optimize import curve_fit

    q = (low - c) / k
    sel = (q > q_min) & np.isfinite(q) & (high > 0)
    if sel.sum() < 2000:
        return {"n_hits": int(sel.sum())}
    q, h = q[sel], high[sel]
    edges = np.linspace(q_min, np.percentile(q, 99.5), 40)
    idx = np.digitize(q, edges) - 1
    xs, ys, ws = [], [], []
    for i in range(len(edges) - 1):
        m = idx == i
        if m.sum() >= 30:
            xs.append(0.5 * (edges[i] + edges[i + 1]))
            ys.append(float(np.median(h[m])))
            ws.append(float(np.sqrt(m.sum())))
    xs, ys, ws = np.array(xs), np.array(ys), np.array(ws)
    if len(xs) < 6:
        return {"n_hits": int(sel.sum())}
    try:
        popt, pcov = curve_fit(roll_over, xs, ys, p0=[2100.0, 8.0],
                               sigma=1.0 / ws, absolute_sigma=False, maxfev=20000,
                               bounds=([500.0, 1.0], [6000.0, 40.0]))
    except Exception as exc:                      # noqa: BLE001
        return {"n_hits": int(sel.sum()), "error": str(exc)}
    a, n = popt
    resid = ys - roll_over(xs, a, n)
    return {"n_hits": int(sel.sum()), "A": float(a), "n": float(n),
            "A_err": float(np.sqrt(pcov[0, 0])), "n_err": float(np.sqrt(pcov[1, 1])),
            "profile": (xs, ys), "rms_resid": float(np.sqrt(np.mean(resid ** 2)))}


def fit_low_gain_scatter(high: np.ndarray, low: np.ndarray, k: float, c: float,
                         hg_min: float = 30.0, hg_max: float = 1400.0) -> dict:
    """Measure how far the data's low gain scatters around the anchor line.

    AdcDigitizer draws the low gain as `k*q + c`, a line; the data put every hit
    within a band whose width has two parts: an additive noise that is the same
    for a 30 ADC hit and a 1000 ADC one (the low-gain shaper's own pedestal
    width), and a part growing with the amplitude (the channel-to-channel
    dispersion of the gain ratio, 5% per channel on run 13, indistinguishable
    from a per-hit spread in a spectrum).  Both are fitted here as
    sigma_lg^2 = noise^2 + (spread * k * hg)^2 on the robust width (half the
    16-84 percentile range, so the 3% of hits on stray gain ratios do not
    count) of `hit_lg - c - k*hit_hg` in bins of hit_hg.  They become
    `lg_noise_adc` / `lg_gain_spread` in digi_calibration.yml.
    """
    sel = (high > hg_min) & (high < hg_max) & np.isfinite(low)
    if sel.sum() < 2000:
        return {"n_hits": int(sel.sum())}
    h, resid = high[sel], (low[sel] - c) - k * high[sel]
    edges = np.array([30, 60, 100, 200, 400, 700, 1000, 1300, 1450], dtype=float)
    idx = np.digitize(h, edges) - 1
    xs, sig, ws = [], [], []
    for i in range(len(edges) - 1):
        m = idx == i
        if m.sum() >= 200:
            q16, q84 = np.percentile(resid[m], [16, 84])
            xs.append(float(np.mean(h[m])))
            sig.append(0.5 * (q84 - q16))
            ws.append(float(m.sum()))
    xs, sig, ws = np.array(xs), np.array(sig), np.array(ws)
    if len(xs) < 3:
        return {"n_hits": int(sel.sum())}
    # linear least squares in sigma^2 = a^2 + b^2 * (k*hg)^2
    design = np.vstack([np.ones_like(xs), (k * xs) ** 2]).T * np.sqrt(ws)[:, None]
    coef, *_ = np.linalg.lstsq(design, sig ** 2 * np.sqrt(ws), rcond=None)
    noise = float(np.sqrt(max(coef[0], 0.0)))
    spread = float(np.sqrt(max(coef[1], 0.0)))
    return {"n_hits": int(sel.sum()), "noise_adc": noise, "gain_spread": spread,
            "profile": (xs, sig)}


def measure_knee(high: np.ndarray, saturation_adc: float) -> dict:
    """Where the high gain stops following the charge, from its own spectrum.

    The ceiling is taken as the **99.99th percentile of all unmasked hits**: the
    point where the population ends.  Measured on the test beam it is ~2050-2150
    ADC, with only single stragglers out to ~2800 -- so the preamp does not hard
    clip at one value, it rolls over, and a single `adc_high_max` is an
    approximation of that roll-over rather than a physical edge.  Using a high
    percentile of the hits ABOVE the saturation threshold instead would land in
    the bulk (most of those hits sit between 1500 and 2000), and the single
    maximum is one hit and one fluctuation, so neither is used.
    """
    positive = high[high > 0]
    if len(positive) < 500:
        return {"n": int(len(positive))}
    above = high[high > saturation_adc]
    return {"n": int(len(above)),
            "p9999": float(np.percentile(positive, 99.99)),
            "p99999": float(np.percentile(positive, 99.999)),
            "max": float(positive.max())}


def _style(ax, xlabel, ylabel, title=None):
    ax.set_xlabel(xlabel, color=C_INK)
    ax.set_ylabel(ylabel, color=C_INK)
    if title:
        ax.set_title(title, color=C_INK, fontsize=11, pad=10)
    ax.grid(True, color=C_GRID, linewidth=0.6)
    ax.set_axisbelow(True)
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)
    for spine in ("left", "bottom"):
        ax.spines[spine].set_color(C_GRID)
    ax.tick_params(colors=C_INK_SOFT)


def plot_lg_spectrum(data, sim, args, outdir):
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(8.2, 5.2))
    hi = float(np.percentile(data[data > 0], 99.9)) if np.any(data > 0) else 300.0
    bins = np.linspace(-20, hi, 180)
    for values, colour, label, n in ((data, C_DATA, "data", args.n_data),
                                     (sim, C_SIM, "simulation", args.n_sim)):
        if len(values) == 0:
            continue
        ax.hist(values, bins=bins, histtype="step", linewidth=1.8, color=colour,
                weights=np.full(len(values), 1.0 / max(n, 1)), label=label)
    ax.set_yscale("log")
    ax.set_xlim(bins[0], bins[-1])
    _style(ax, "hit_lg  [ADC, low gain, pedestal-subtracted]",
           "hits per event / bin", f"Low-gain spectrum — {args.title}")
    ax.legend(frameon=False, labelcolor=C_INK)
    fig.tight_layout()
    path = os.path.join(outdir, f"lg_spectrum_{args.tag}.png")
    fig.savefig(path, dpi=140)
    plt.close(fig)
    return path


def plot_lg_vs_hg(data_hg, data_lg, sim_hg, sim_lg, args, outdir, fit=None):
    import matplotlib.pyplot as plt

    anchor = ANCHOR.get(args.threshold, ANCHOR["th230"])
    hi_x = float(np.percentile(data_hg, 99.99)) if len(data_hg) else 3000.0
    hi_x = max(hi_x, args.saturation_adc * 2)
    hi_y = anchor["k"] * hi_x + anchor["c"]

    fig, axes = plt.subplots(1, 2, figsize=(11.6, 5.0), sharex=True, sharey=True)
    for ax, (hg, lg, label) in zip(axes, ((data_hg, data_lg, "data"),
                                          (sim_hg, sim_lg, "simulation"))):
        sel = (hg > 0) & (lg != 0)
        if sel.sum():
            ax.hexbin(hg[sel], lg[sel], gridsize=90, bins="log", cmap="BuPu",
                      mincnt=1, extent=(0, hi_x, 0, hi_y))
        line = np.linspace(0, hi_x, 200)
        ax.plot(line, anchor["k"] * line + anchor["c"], color=C_INK_SOFT,
                linewidth=1.4, linestyle="--",
                label=f"anchor k={anchor['k']:.4f}, c={anchor['c']:.2f}")
        ax.axvline(args.saturation_adc, color=C_SIM, linewidth=1.1, linestyle=":")
        if fit and "A" in fit:
            # The fitted roll-over, drawn in the LG-vs-HG plane: a charge q gives
            # hit_lg = k*q + c and hit_hg = S(q).
            qq = np.linspace(0, hi_x * 1.5, 400)
            ax.plot(roll_over(qq, fit["A"], fit["n"]), anchor["k"] * qq + anchor["c"],
                    color=C_SIM, linewidth=1.6,
                    label=f"roll-over fit A={fit['A']:.0f}, n={fit['n']:.1f}")
        _style(ax, "hit_hg  [ADC]", "hit_lg  [ADC]", label)
        ax.set_xlim(0, hi_x)
        ax.set_ylim(0, hi_y)
    axes[0].legend(frameon=False, labelcolor=C_INK, loc="upper left", fontsize=9)
    fig.suptitle(f"Low gain against high gain — {args.title}", color=C_INK,
                 fontsize=11)
    fig.tight_layout()
    path = os.path.join(outdir, f"lg_vs_hg_{args.tag}.png")
    fig.savefig(path, dpi=140)
    plt.close(fig)
    return path


def main(argv=None) -> int:
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--data", required=True, help="Test-beam ecal tree")
    p.add_argument("--sim", required=True,
                   help="Simulated ecal tree made WITH the ADC model")
    p.add_argument("--threshold", default="th230", choices=sorted(ANCHOR),
                   help="Threshold set, for the anchor line drawn on top")
    p.add_argument("--saturation-adc", type=float, default=1500.0)
    p.add_argument("--max-data-events", type=int, default=20000)
    p.add_argument("--max-sim-events", type=int, default=None)
    p.add_argument("--write", action="store_true",
                   help="Write the fitted A / n into mappings/digi_calibration.yml "
                        "as adc_high_max / saturation_order for this threshold set")
    p.add_argument("--tag", default="gains")
    p.add_argument("--title", default="")
    p.add_argument("--outdir", default="plots")
    args = p.parse_args(argv)
    if not args.title:
        args.title = args.tag

    import matplotlib
    matplotlib.use("Agg")

    os.makedirs(args.outdir, exist_ok=True)
    data_hg, data_lg, args.n_data = read_hits(args.data, args.max_data_events)
    sim_hg, sim_lg, args.n_sim = read_hits(args.sim, args.max_sim_events)

    knee_data = measure_knee(data_hg, args.saturation_adc)
    knee_sim = measure_knee(sim_hg, args.saturation_adc)
    anchor = ANCHOR[args.threshold]
    fit = fit_dynamic_range(data_hg, data_lg, anchor["k"], anchor["c"])
    scatter = fit_low_gain_scatter(data_hg, data_lg, anchor["k"], anchor["c"])
    scatter_sim = fit_low_gain_scatter(sim_hg, sim_lg, anchor["k"], anchor["c"])

    written = [plot_lg_spectrum(data_lg, sim_lg, args, args.outdir),
               plot_lg_vs_hg(data_hg, data_lg, sim_hg, sim_lg, args, args.outdir,
                             fit=fit)]

    lines = [f"# Low gain vs high gain — {args.title}",
             f"#   threshold set {args.threshold}, anchor "
             f"k={ANCHOR[args.threshold]['k']}, c={ANCHOR[args.threshold]['c']} "
             f"(from {ANCHOR[args.threshold]['source_run']})",
             "",
             f"{'':<14}{'hits > sat':>12}{'p99.99':>9}{'p99.999':>9}{'max':>9}"]
    for label, knee in (("data", knee_data), ("simulation", knee_sim)):
        if "p9999" in knee:
            lines.append(f"{label:<14}{knee['n']:>12d}{knee['p9999']:>9.0f}"
                         f"{knee['p99999']:>9.0f}{knee['max']:>9.0f}")
        else:
            lines.append(f"{label:<14}{knee.get('n', 0):>12d}"
                         f"{'too few':>9}{'':>9}{'':>9}")
    if "A" in fit:
        lines += ["",
                  f"dynamic range fitted from the data's LG-vs-HG relation "
                  f"({fit['n_hits']} hits above q_min):",
                  f"  S(q) = q / (1 + (q/A)^n)^(1/n)   A = {fit['A']:.0f} +- "
                  f"{fit['A_err']:.0f} ADC   n = {fit['n']:.2f} +- {fit['n_err']:.2f}"
                  f"   rms residual {fit['rms_resid']:.1f} ADC",
                  f"  (population edge, p99.99 of hit_hg: {knee_data.get('p9999', float('nan')):.0f} ADC)"]
        if "noise_adc" in scatter:
            lines += ["",
                      "low-gain scatter around the anchor line, robust width of "
                      "hit_lg - c - k*hit_hg in bins of hit_hg (30-1400 ADC):",
                      f"  data       noise {scatter['noise_adc']:.2f} LG-ADC, gain spread "
                      f"{100 * scatter['gain_spread']:.1f} %  ({scatter['n_hits']} hits)"]
            if "noise_adc" in scatter_sim:
                lines.append(f"  simulation noise {scatter_sim['noise_adc']:.2f} LG-ADC, "
                             f"gain spread {100 * scatter_sim['gain_spread']:.1f} %")
        if args.write:
            import yaml
            calib_file = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                      "..", "mappings", "digi_calibration.yml")
            with open(calib_file) as fh:
                calib = yaml.safe_load(fh)
            entry = calib["thresholds"][args.threshold]
            entry["adc_high_max"] = round(fit["A"], 1)
            entry["saturation_order"] = round(fit["n"], 2)
            entry["adc_high_max_source"] = (
                f"roll-over fitted to {os.path.basename(args.data)} "
                f"(hit_lg-vs-hit_hg, {fit['n_hits']} hits)")
            entry["gain_ratio"] = anchor["k"]
            entry["gain_intercept"] = anchor["c"]
            if "noise_adc" in scatter:
                entry["lg_noise_adc"] = round(scatter["noise_adc"], 2)
                entry["lg_gain_spread"] = round(scatter["gain_spread"], 4)
                entry["lg_scatter_source"] = (
                    f"robust width of hit_lg about the anchor line in bins of hit_hg, "
                    f"{os.path.basename(args.data)} ({scatter['n_hits']} hits, 30-1400 ADC)")
            with open(calib_file, "w") as fh:
                yaml.safe_dump(calib, fh, sort_keys=False, allow_unicode=True)
            lines.append(f"  written to {os.path.normpath(calib_file)}")
    elif "error" in fit:
        lines += ["", f"dynamic-range fit failed: {fit['error']}"]
    summary = "\n".join(lines)
    print("\n" + summary)
    summary_path = os.path.join(args.outdir, f"summary_gains_{args.tag}.txt")
    with open(summary_path, "w") as fh:
        fh.write(summary + "\n")
    written.append(summary_path)

    print("\nWritten:")
    for path in written:
        print(f"  {path}")
    return 0


if __name__ == "__main__":
    sys.exit(main())

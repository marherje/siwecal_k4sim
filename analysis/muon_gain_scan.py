#!/usr/bin/env python3
"""
Gain scan on a muon run: which `adc_per_mip` makes the simulated per-hit
high-gain spectrum look like the data's, using the whole spectrum shape.

Why a scan and not the fixed-point fit of `fit_adc_scale.py --observable
hit-peak`: at th230 the discriminator (25.4 ADC) sits ABOVE the MIP, so what
survives is the Landau tail, which is scale-free -- the mode of the surviving
spectrum is not monotonic in the gain (26.6 / 27.6 / 27.5 / 28.4 / 26.6 ADC for
gain 20 / 21 / 22 / 23 / 24) and a fixed-point iteration on it lands wherever
it started.  The shape just above threshold still knows where the MPV is: the
spectrum falls from the threshold when the MPV is below it and rises first when
it is above.  A chi2 of the unit-normalised spectra over [lo, hi] ADC against
the gain has one minimum, and its position is the gain.

The data's own event selection is applied to both sides: `EcalEventBuilder`
keeps an event only if `MinSlabsHit = 10` slabs fired (a Landau selection on a
muon's fifteen crossings), and events above `--max-hits` are pile-up.

Usage
-----
    python -m analysis.muon_gain_scan --threshold th230 \\
        --data /eos/.../ecal_TB2026CERN_run_000004.root \\
        --sim 20=mu_g20.root --sim 21=mu_g21.root ... --outdir plots/
"""

from __future__ import annotations

import argparse
import os
import sys

import numpy as np

from analysis.fit_adc_scale import hit_mip_peak, layers_hit


def selected_spectrum(path, lo, hi, max_hits, min_layers, max_events=None):
    import uproot

    tree = uproot.open(path)["ecal"]
    arrays = tree.arrays(["nhit_chan", "hit_slab", "hit_hg", "hit_ismasked"],
                         entry_stop=max_events, library="np")
    keep = arrays["nhit_chan"] <= max_hits
    keep &= layers_hit(arrays["hit_slab"], arrays["hit_ismasked"]) >= min_layers
    high = np.concatenate(arrays["hit_hg"][keep]).astype(float)
    masked = np.concatenate(arrays["hit_ismasked"][keep]).astype(int)
    high = high[masked == 0]
    counts, edges = np.histogram(high, bins=np.arange(0, 151, 1))
    norm = counts[lo:hi].sum()
    return counts / max(norm, 1), norm, int(keep.sum())


def main(argv=None) -> int:
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--threshold", required=True)
    p.add_argument("--data", required=True, help="Muon run ecal tree")
    p.add_argument("--sim", action="append", default=[], metavar="GAIN=PATH",
                   help="Digitised muon sample at that adc_per_mip, repeatable")
    p.add_argument("--lo", type=int, default=None,
                   help="Lower edge of the chi2 window [ADC] "
                        "(default: threshold_adc + 1)")
    p.add_argument("--hi", type=int, default=80)
    p.add_argument("--max-hits", type=int, default=30)
    p.add_argument("--min-layers", type=int, default=10)
    p.add_argument("--max-data-events", type=int, default=40000)
    p.add_argument("--outdir", default="plots")
    p.add_argument("--tag", default="")
    args = p.parse_args(argv)

    import yaml
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    calib_file = os.path.join(os.path.dirname(__file__), "..", "mappings",
                              "digi_calibration.yml")
    with open(calib_file) as fh:
        entry = yaml.safe_load(fh)["thresholds"][args.threshold]
    lo = args.lo if args.lo is not None else int(round(entry["threshold_adc"])) + 1
    hi = args.hi
    tag = args.tag or args.threshold
    os.makedirs(args.outdir, exist_ok=True)

    d_spec, d_norm, d_events = selected_spectrum(
        args.data, lo, hi, args.max_hits, args.min_layers, args.max_data_events)
    d_peak = hit_mip_peak(args.data, args.max_data_events,
                          max_hits=args.max_hits, min_layers=args.min_layers)
    print(f"[data] {os.path.basename(args.data)}: {d_events} events selected "
          f"(<= {args.max_hits} hits, >= {args.min_layers} slabs), "
          f"{d_norm} hits in [{lo}, {hi}) ADC, peak {d_peak:.1f} ADC")

    rows = []
    for spec in args.sim:
        gain, _, path = spec.partition("=")
        gain = float(gain)
        s_spec, s_norm, s_events = selected_spectrum(
            path, lo, hi, args.max_hits, args.min_layers)
        # Pearson chi2 of the unit-normalised spectra, data as the expectation.
        chi2 = float(np.sum((s_spec[lo:hi] - d_spec[lo:hi]) ** 2
                            / np.maximum(d_spec[lo:hi], 1e-4)))
        peak = hit_mip_peak(path, max_hits=args.max_hits, min_layers=args.min_layers)
        rows.append((gain, chi2, peak, s_events, s_norm, s_spec))
        print(f"[sim ] gain {gain:5.1f}: {s_events:5d} events selected, "
              f"peak {peak:.1f} ADC, chi2 {chi2 * 1e3:.1f}e-3")
    rows.sort()
    gains = np.array([r[0] for r in rows])
    chi2s = np.array([r[1] for r in rows])

    # Parabola through the minimum and its neighbours: the gain, and the +-1
    # (in units of the minimum chi2) half-width as the interval.
    best = None
    if len(rows) >= 3:
        i = int(np.argmin(chi2s))
        i = min(max(i, 1), len(rows) - 2)
        a, b, c = np.polyfit(gains[i - 1:i + 2], chi2s[i - 1:i + 2], 2)
        if a > 0:
            g0 = -b / (2 * a)
            chi2_min = c - b * b / (4 * a)
            half = np.sqrt(chi2_min / a) if chi2_min > 0 else float("nan")
            best = (g0, half)
            print(f"[scan] minimum at adc_per_mip = {g0:.2f} (+-{half:.2f} where "
                  f"chi2 doubles)")

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 4.6))
    x = np.arange(0, 150) + 0.5
    ax1.step(x, d_spec, where="mid", color="black", linewidth=2.0,
             label=f"data ({os.path.basename(args.data).split('.')[0][-6:]})")
    cmap = plt.get_cmap("viridis")
    for k, (gain, chi2, peak, s_events, s_norm, s_spec) in enumerate(rows):
        ax1.step(x, s_spec, where="mid", color=cmap(k / max(len(rows) - 1, 1)),
                 linewidth=1.2, label=f"sim (digi), {gain:g} ADC/MIP")
    ax1.axvspan(lo, hi, color="grey", alpha=0.10, label=f"$\\chi^2$ window [{lo}, {hi})")
    ax1.axvline(entry["threshold_adc"], color="red", linestyle=":", linewidth=1.0,
                label=f"discriminator {entry['threshold_adc']:.1f} ADC")
    ax1.set_xlim(0, 120)
    ax1.set_xlabel("hit high gain, pedestal subtracted [ADC]")
    ax1.set_ylabel(f"fraction of hits in [{lo}, {hi}) per 1 ADC")
    ax1.set_title(f"{args.threshold} muons — spectrum shape vs gain")
    ax1.legend(fontsize=7.5)
    ax1.grid(alpha=0.3)

    ax2.plot(gains, chi2s * 1e3, "o-", color="tab:blue")
    if best is not None:
        ax2.axvline(best[0], color="tab:red", linestyle="--",
                    label=f"minimum: {best[0]:.2f} ± {best[1]:.2f} ADC/MIP")
        ax2.legend()
    ax2.set_xlabel("adc_per_mip [ADC/MIP]")
    ax2.set_ylabel("$\\chi^2$ of the normalised spectra ×10³")
    ax2.set_title("shape $\\chi^2$ against the gain")
    ax2.grid(alpha=0.3)
    fig.tight_layout()
    out = os.path.join(args.outdir, f"muon_gain_scan_{tag}.png")
    fig.savefig(out, dpi=130)
    print(f"Written: {out}")

    with open(os.path.join(args.outdir, f"muon_gain_scan_{tag}.txt"), "w") as fh:
        fh.write(f"# {args.threshold}: chi2 window [{lo}, {hi}) ADC, events <= "
                 f"{args.max_hits} hits and >= {args.min_layers} slabs\n")
        fh.write(f"# data peak {d_peak:.1f} ADC, {d_events} events\n")
        fh.write("gain  chi2x1e3  sim_peak  sim_events\n")
        for gain, chi2, peak, s_events, s_norm, _ in rows:
            fh.write(f"{gain:5.1f} {chi2 * 1e3:9.2f} {peak:9.1f} {s_events:11d}\n")
        if best is not None:
            fh.write(f"minimum {best[0]:.3f} +- {best[1]:.3f}\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())

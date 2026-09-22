#!/usr/bin/env python3
"""
Compare the ADC response of the digitised simulation with test-beam data.

Why a conversion is needed
--------------------------
The two ecal trees do NOT carry the same quantity:

  data (siwecal-tb2026, EcalEventBuilder)
      hit_hg      = adc_high - pedestal_hg   [ADC]  <-- the raw detector reading
      hit_energy  = hit_hg / MIP_hg          [MIP]  (low-gain-anchored above
                                                     AdcSaturationThreshold)
  simulation (analysis/sim_to_ecal_tree.py)
      hit_hg      = 0                               <-- the simulation has no ADC
      hit_energy  = digitised amplitude      [MIP]

so the simulated MIP are converted back to ADC with the SAME per-channel MIP
table the data reconstruction used (mips/th<N>/MIP_*_highgain.txt):

    ADC_sim(hit) = hit_energy[MIP] * mpv(slab, chip, channel)

The channel indices exist in the simulation because ChannelMapper (job 3) has
already rewritten the CellIDs into the test-beam format, and the masking flags
come from that very same table, so both samples are missing the same channels.

Saturation
----------
Above ``--saturation-adc`` (1500, the event builder's AdcSaturationThreshold)
the data's hit_hg is flat -- the high-gain preamp is saturated -- while its
hit_energy is recovered from the low-gain branch. The simulation has NO
saturation model, so two data curves are drawn:

    "data (hit_hg)"       the raw pedestal-subtracted ADC, saturation included
    "data (linearised)"   hit_energy * mpv, i.e. the ADC an unsaturated preamp
                          would have given -- the apples-to-apples comparison

Usage
-----
    python -m analysis.compare_adc_data_sim \\
        --data  /eos/.../Reconstruction/TB2026CERN_run_000013/ecal_TB2026CERN_run_000013.root \\
        --sim   digi=/path/ecal_sim_e52_real.root \\
        --sim   simple=/path/ecal_sim_e52_simple.root \\
        --mip-file ../siwecal-tb2026/calibration/MuonCalib_gaudi/mips/th230/\\
MIP_pedestalsubmode1_TB2026CERN_run_000004_highgain.txt \\
        --tag e52 --title "e- 52 GeV  vs  TB2026CERN run 13 (P1, 52 GeV)" \\
        --outdir plots/
"""

from __future__ import annotations

import argparse
import os
import sys
from dataclasses import dataclass, field
from typing import Dict, List, Optional

import numpy as np

# --------------------------------------------------------------------------- #
# Detector constants (kept in step with siwecal-tb2026)
# --------------------------------------------------------------------------- #
N_LAYERS = 15
N_CHIPS = 16
N_CHANS = 64

# ChannelMapper.MaxMIPValue: a channel whose MPV is 0, NaN or above this is not
# calibrated, and is masked in both samples.
MAX_MIP_ADC = 100.0

# --------------------------------------------------------------------------- #
# Palette: the data is ONE entity read two ways (solid / dashed, same hue); each
# digitisation chain gets its own hue, assigned in fixed order.
# --------------------------------------------------------------------------- #
C_DATA = "#14407f"
C_SIM = ["#b03024", "#e8856b", "#eda100", "#e87ba4"]
C_GRID = "#d8d8d4"
C_INK = "#0b0b0b"
C_INK_SOFT = "#52514e"


# --------------------------------------------------------------------------- #
# MIP calibration table
# --------------------------------------------------------------------------- #

def load_mip_table(path: str) -> np.ndarray:
    """(layer, chip, channel) -> MPV [ADC/MIP]; NaN where the channel is masked.

    Same file and same masking rule as ChannelMapper (mpv == 0, NaN, or > 100).
    """
    raw = np.loadtxt(path, comments="#")
    if raw.ndim != 2 or raw.shape[1] < 4:
        sys.exit(f"ERROR: {path} is not a 'layer chip channel mpv ...' table")

    table = np.full((N_LAYERS, N_CHIPS, N_CHANS), np.nan)
    layer = raw[:, 0].astype(int)
    chip = raw[:, 1].astype(int)
    chan = raw[:, 2].astype(int)
    mpv = raw[:, 3].astype(float)

    ok = (mpv > 0.0) & (mpv <= MAX_MIP_ADC) & np.isfinite(mpv)
    ok &= (layer >= 0) & (layer < N_LAYERS)
    ok &= (chip >= 0) & (chip < N_CHIPS)
    ok &= (chan >= 0) & (chan < N_CHANS)
    table[layer[ok], chip[ok], chan[ok]] = mpv[ok]

    n_good = int(np.isfinite(table).sum())
    print(f"[mip] {os.path.basename(path)}: {n_good} calibrated channels "
          f"({n_good / (N_LAYERS * N_CHIPS * N_CHANS) * 100:.1f}%), "
          f"MPV mean {np.nanmean(table):.2f} ADC/MIP")
    return table


def mpv_of(table: np.ndarray, slab, chip, chan) -> np.ndarray:
    """Per-hit MPV lookup; NaN for out-of-range or uncalibrated channels."""
    ok = ((slab >= 0) & (slab < N_LAYERS) & (chip >= 0) & (chip < N_CHIPS)
          & (chan >= 0) & (chan < N_CHANS))
    out = np.full(len(slab), np.nan)
    out[ok] = table[slab[ok], chip[ok], chan[ok]]
    return out


# --------------------------------------------------------------------------- #
# Accumulated histograms (the data run has O(100k) events: never load it whole)
# --------------------------------------------------------------------------- #

@dataclass
class Sample:
    """Histograms and running sums for one sample."""

    label: str
    color: str
    dashed: bool = False
    n_events: int = 0
    n_hits: int = 0
    h_hit_adc: np.ndarray = field(default=None)
    h_sum_adc: np.ndarray = field(default=None)
    h_nhit: np.ndarray = field(default=None)
    layer_adc: np.ndarray = field(default_factory=lambda: np.zeros(N_LAYERS))
    layer_hits: np.ndarray = field(default_factory=lambda: np.zeros(N_LAYERS))
    sum_adc_values: List[np.ndarray] = field(default_factory=list)
    hit_adc_sum: float = 0.0
    n_saturated: int = 0

    def init_hists(self, bins):
        self.h_hit_adc = np.zeros(len(bins.hit) - 1)
        self.h_sum_adc = np.zeros(len(bins.sum) - 1)
        self.h_nhit = np.zeros(len(bins.nhit) - 1)


@dataclass
class Bins:
    hit: np.ndarray
    sum: np.ndarray
    nhit: np.ndarray


def fill(sample: Sample, bins: Bins, adc, slab, keep, counts,
         saturation_adc: float):
    """Add one chunk.

    ``adc``/``slab``/``keep`` are flat per-hit arrays and ``counts`` the
    per-event hit multiplicities that slice them.  The event index is built with
    ``np.repeat``, so events left with no hit by the cut still count as events
    (with sum 0) instead of shifting everything that follows.
    """
    sample.n_events += len(counts)
    if len(counts) == 0:
        return

    event_idx = np.repeat(np.arange(len(counts)), counts)
    kept_adc = adc[keep]
    kept_slab = slab[keep]

    sample.n_hits += int(keep.sum())
    sample.hit_adc_sum += float(kept_adc.sum())
    sample.n_saturated += int((kept_adc >= saturation_adc).sum())
    sample.h_hit_adc += np.histogram(kept_adc, bins=bins.hit)[0]

    in_range = (kept_slab >= 0) & (kept_slab < N_LAYERS)
    sample.layer_adc += np.bincount(kept_slab[in_range],
                                    weights=kept_adc[in_range],
                                    minlength=N_LAYERS)[:N_LAYERS]
    sample.layer_hits += np.bincount(kept_slab[in_range],
                                     minlength=N_LAYERS)[:N_LAYERS]

    sums = np.bincount(event_idx, weights=np.where(keep, adc, 0.0),
                       minlength=len(counts))
    n_kept = np.bincount(event_idx, weights=keep.astype(np.float64),
                         minlength=len(counts))

    sample.h_sum_adc += np.histogram(sums, bins=bins.sum)[0]
    sample.h_nhit += np.histogram(n_kept, bins=bins.nhit)[0]
    sample.sum_adc_values.append(sums.astype(np.float32))


# --------------------------------------------------------------------------- #
# Readers
# --------------------------------------------------------------------------- #

_BRANCHES = ["nhit_chan", "hit_slab", "hit_chip", "hit_chan",
             "hit_energy", "hit_ismasked"]


def read_sample(path: str, table: np.ndarray, bins: Bins, label: str,
                color: str, *, use_hit_hg: bool, mip_cut: float,
                saturation_adc: float, max_events: Optional[int],
                tree_name: str = "ecal", dashed: bool = False,
                step: int = 20000) -> Sample:
    """Stream one ecal tree, filling ``Sample``.

    ``use_hit_hg`` picks the RAW data ADC (hit_hg, saturation included); other-
    wise the ADC is hit_energy * MPV -- the linearised data reading, and the only
    way to get an ADC out of the simulation at all.
    """
    import uproot

    sample = Sample(label=label, color=color, dashed=dashed)
    sample.init_hists(bins)

    branches = list(_BRANCHES) + (["hit_hg"] if use_hit_hg else [])
    tree = uproot.open(path)[tree_name]
    n_total = tree.num_entries
    stop = min(n_total, max_events) if max_events else n_total

    def flat(chunk, name, dtype):
        parts = chunk[name]
        if len(parts) == 0:
            return np.zeros(0, dtype=dtype)
        return np.concatenate(parts).astype(dtype)

    for chunk in tree.iterate(branches, entry_stop=stop, step_size=step,
                              library="np"):
        counts = chunk["nhit_chan"].astype(np.int64)
        slab = flat(chunk, "hit_slab", np.int64)
        chip = flat(chunk, "hit_chip", np.int64)
        chan = flat(chunk, "hit_chan", np.int64)
        energy = flat(chunk, "hit_energy", np.float64)
        masked = flat(chunk, "hit_ismasked", np.int64)

        # nhit_chan is what slices the hit arrays; a mismatch would silently
        # pair hits with the wrong events.
        if counts.sum() != len(slab):
            sys.exit(f"ERROR: {path}: nhit_chan sums to {counts.sum()} but "
                     f"hit_slab has {len(slab)} entries")

        mpv = mpv_of(table, slab, chip, chan)
        adc = flat(chunk, "hit_hg", np.float64) if use_hit_hg else energy * mpv

        # One selection for every sample: calibrated channel, same MIP cut.
        # np.isfinite(mpv) is required even on the raw-ADC curve, so that the
        # two data readings are drawn from exactly the same set of hits.
        keep = ((masked == 0) & np.isfinite(mpv) & np.isfinite(adc)
                & (energy > mip_cut))

        fill(sample, bins, adc, slab, keep, counts, saturation_adc)

    print(f"[read] {label:<22} {sample.n_events:7d} events, "
          f"{sample.n_hits:9d} hits  ({os.path.basename(path)})")
    return sample


# --------------------------------------------------------------------------- #
# Plots
# --------------------------------------------------------------------------- #

def _style(ax, xlabel, ylabel, title=None):
    ax.set_xlabel(xlabel, color=C_INK)
    ax.set_ylabel(ylabel, color=C_INK)
    if title:
        ax.set_title(title, color=C_INK, fontsize=11, pad=10)
    ax.grid(True, color=C_GRID, linewidth=0.6, alpha=0.9)
    ax.set_axisbelow(True)
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)
    for spine in ("left", "bottom"):
        ax.spines[spine].set_color(C_GRID)
    ax.tick_params(colors=C_INK_SOFT)


def _step(ax, edges, values, sample: Sample, scale=1.0):
    # stairs() closes the last bin; ax.step(edges[:-1], ...) would drop it.
    ax.stairs(values * scale, edges, color=sample.color, linewidth=1.6,
              linestyle="--" if sample.dashed else "-", label=sample.label)


def plot_hit_spectrum(samples, bins, args, outdir):
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(8.2, 5.2))
    for s in samples:
        _step(ax, bins.hit, s.h_hit_adc / max(s.n_events, 1), s)
    ax.set_xlim(bins.hit[0], bins.hit[-1])
    ax.set_yscale("log")
    # Only worth drawing when the axis actually reaches the saturation point.
    if bins.hit[0] < args.saturation_adc < bins.hit[-1]:
        ax.axvline(args.saturation_adc, color=C_INK_SOFT, linewidth=1.2,
                   linestyle=":", zorder=0)
        ax.annotate("high gain saturates\n(AdcSaturationThreshold)",
                    xy=(args.saturation_adc, ax.get_ylim()[1]),
                    xytext=(-8, -12), textcoords="offset points",
                    ha="right", va="top", fontsize=8, color=C_INK_SOFT)
    _style(ax, "hit amplitude  [ADC, high gain, pedestal-subtracted]",
           "hits per event / bin",
           f"Per-hit ADC spectrum — {args.title}")
    ax.legend(frameon=False, labelcolor=C_INK)
    fig.tight_layout()
    path = os.path.join(outdir, f"hit_adc_spectrum_{args.tag}.png")
    fig.savefig(path, dpi=140)
    plt.close(fig)
    return path


def plot_event_sum(samples, bins, args, outdir):
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(8.2, 5.2))
    for s in samples:
        area = s.h_sum_adc.sum() * np.diff(bins.sum)[0]
        if area > 0:
            _step(ax, bins.sum, s.h_sum_adc / area, s)
    ax.set_xlim(bins.sum[0], bins.sum[-1])
    _style(ax, "event total  [ADC]", "fraction of events / ADC",
           f"Event ADC sum — {args.title}")
    ax.legend(frameon=False, labelcolor=C_INK)
    fig.tight_layout()
    path = os.path.join(outdir, f"event_sum_adc_{args.tag}.png")
    fig.savefig(path, dpi=140)
    plt.close(fig)
    return path


def plot_layer_profile(samples, args, outdir):
    import matplotlib.pyplot as plt

    layers = np.arange(N_LAYERS)
    fig, (ax, axh) = plt.subplots(2, 1, figsize=(8.2, 7.0), sharex=True,
                                  gridspec_kw={"height_ratios": [2, 1]})
    for s in samples:
        n = max(s.n_events, 1)
        ax.plot(layers, s.layer_adc / n, color=s.color, linewidth=2.0,
                linestyle="--" if s.dashed else "-", marker="o",
                markersize=4, label=s.label)
        axh.plot(layers, s.layer_hits / n, color=s.color, linewidth=2.0,
                 linestyle="--" if s.dashed else "-", marker="o", markersize=4)
    _style(ax, "", "ADC per event", f"Longitudinal profile — {args.title}")
    _style(axh, "layer", "hits per event")
    ax.legend(frameon=False, labelcolor=C_INK)
    fig.tight_layout()
    path = os.path.join(outdir, f"layer_profile_{args.tag}.png")
    fig.savefig(path, dpi=140)
    plt.close(fig)
    return path


def plot_nhit(samples, bins, args, outdir):
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(8.2, 5.2))
    for s in samples:
        area = s.h_nhit.sum() * np.diff(bins.nhit)[0]
        if area > 0:
            _step(ax, bins.nhit, s.h_nhit / area, s)
    ax.set_xlim(bins.nhit[0], bins.nhit[-1])
    _style(ax, "hits per event (above the common MIP cut)",
           "fraction of events / hit",
           f"Hit multiplicity — {args.title}")
    ax.legend(frameon=False, labelcolor=C_INK)
    fig.tight_layout()
    path = os.path.join(outdir, f"nhit_{args.tag}.png")
    fig.savefig(path, dpi=140)
    plt.close(fig)
    return path


# --------------------------------------------------------------------------- #
# Summary
# --------------------------------------------------------------------------- #

def write_layer_csv(samples: List[Sample], args, outdir) -> str:
    """Per-layer ADC and hits per event, plus each sample's ratio to the data."""
    ref = samples[0]
    ref_layer = np.where(ref.layer_adc > 0, ref.layer_adc, np.nan) / max(ref.n_events, 1)

    header = ["layer"]
    for s in samples:
        tag = s.label.replace(" ", "_")
        header += [f"adc_per_event[{tag}]", f"hits_per_event[{tag}]",
                   f"ratio_to_data[{tag}]"]
    rows = []
    for layer in range(N_LAYERS):
        row = [str(layer)]
        for s in samples:
            n = max(s.n_events, 1)
            adc = s.layer_adc[layer] / n
            row += [f"{adc:.1f}", f"{s.layer_hits[layer] / n:.2f}",
                    f"{adc / ref_layer[layer]:.3f}" if ref_layer[layer] > 0 else ""]
        rows.append(",".join(row))

    path = os.path.join(outdir, f"layer_profile_{args.tag}.csv")
    with open(path, "w") as fh:
        fh.write(",".join(header) + "\n")
        fh.write("\n".join(rows) + "\n")
    return path


def summarise(samples: List[Sample], bins: Bins, args) -> str:
    ref = samples[0]
    ref_sum = np.concatenate(ref.sum_adc_values) if ref.sum_adc_values else np.zeros(1)
    ref_mean = float(ref_sum.mean()) if len(ref_sum) else 0.0

    lines = [
        f"# ADC comparison — {args.title}",
        f"#   MIP table   : {args.mip_file}",
        f"#   hit cut     : hit_energy > {args.mip_cut} MIP, unmasked, calibrated channel",
        f"#   saturation  : {args.saturation_adc} ADC (data hit_hg is flat above it)",
        "",
        f"{'sample':<24}{'events':>8}{'hits/ev':>9}{'<ADC>/hit':>11}"
        f"{'peak ADC':>10}{'<sum ADC>':>12}{'sigma/mu':>10}{'/data':>8}{'>sat':>8}",
    ]
    for s in samples:
        n = max(s.n_events, 1)
        sums = np.concatenate(s.sum_adc_values) if s.sum_adc_values else np.zeros(1)
        mu = float(sums.mean())
        sigma = float(sums.std())
        centres = 0.5 * (bins.hit[:-1] + bins.hit[1:])
        peak = float(centres[int(np.argmax(s.h_hit_adc))]) if s.h_hit_adc.sum() else 0.0
        lines.append(
            f"{s.label:<24}{s.n_events:>8d}{s.n_hits / n:>9.1f}"
            f"{s.hit_adc_sum / max(s.n_hits, 1):>11.1f}{peak:>10.1f}"
            f"{mu:>12.0f}{(sigma / mu if mu else 0):>10.3f}"
            f"{(mu / ref_mean if ref_mean else 0):>8.2f}"
            f"{s.n_saturated / max(s.n_hits, 1) * 100:>7.2f}%")
    return "\n".join(lines)


# --------------------------------------------------------------------------- #
# CLI
# --------------------------------------------------------------------------- #

def _parse_args(argv=None):
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--data", required=True, help="Test-beam ecal tree")
    p.add_argument("--sim", action="append", default=[], metavar="NAME=PATH",
                   help="Simulation ecal tree, repeatable (e.g. digi=/path.root). "
                        "NAME is what the legend says, so name the chain, not "
                        "the collection: 'digi' for the RealDigitizer chain, "
                        "'simple' for the threshold-cut one")
    p.add_argument("--mip-file", required=True,
                   help="MIP table used by the data reconstruction "
                        "(mips/th<N>/MIP_*_highgain.txt)")
    p.add_argument("--tag", default="cmp", help="File-name tag for the outputs")
    p.add_argument("--title", default="", help="Plot title")
    p.add_argument("--outdir", default="plots", help="Output directory")
    p.add_argument("--mip-cut", type=float, default=0.5,
                   help="Common per-hit cut [MIP] (default: 0.5, the "
                        "RealDigitizer trigger threshold)")
    p.add_argument("--saturation-adc", type=float, default=1500.0,
                   help="Event builder's AdcSaturationThreshold (default: 1500)")
    p.add_argument("--max-data-events", type=int, default=20000,
                   help="Data events to read (0 = all; default: 20000)")
    p.add_argument("--hit-adc-max", type=float, default=3000.0)
    p.add_argument("--hit-adc-bins", type=int, default=300)
    p.add_argument("--sum-adc-max", type=float, default=None,
                   help="Upper edge of the event-sum axis (default: auto)")
    p.add_argument("--nhit-max", type=float, default=None)
    return p.parse_args(argv)


def main(argv=None) -> int:
    args = _parse_args(argv)
    if not args.sim:
        sys.exit("ERROR: at least one --sim NAME=PATH is required")
    if not args.title:
        args.title = args.tag

    import matplotlib
    matplotlib.use("Agg")

    os.makedirs(args.outdir, exist_ok=True)
    table = load_mip_table(args.mip_file)

    sim_paths: Dict[str, str] = {}
    for spec in args.sim:
        if "=" not in spec:
            sys.exit(f"ERROR: --sim expects NAME=PATH, got '{spec}'")
        name, path = spec.split("=", 1)
        sim_paths[name] = path

    # A first pass over the simulation fixes the axes: the data run is the one
    # that must not be read twice.
    bins = Bins(hit=np.linspace(-50.0, args.hit_adc_max, args.hit_adc_bins + 1),
                sum=np.linspace(0.0, 1.0, 2), nhit=np.linspace(0.0, 1.0, 2))

    import uproot
    sum_guess, nhit_guess = [], []
    for path in sim_paths.values():
        t = uproot.open(path)["ecal"]
        arr = t.arrays(["sum_energy", "nhit_chan"], library="np")
        sum_guess.append(np.percentile(arr["sum_energy"], 99) * float(np.nanmean(table)))
        nhit_guess.append(np.percentile(arr["nhit_chan"], 99))
    sum_max = args.sum_adc_max or float(max(sum_guess)) * 1.6
    nhit_max = args.nhit_max or float(max(nhit_guess)) * 1.6
    bins.sum = np.linspace(0.0, sum_max, 121)
    bins.nhit = np.linspace(0.0, nhit_max, 121)

    max_data = args.max_data_events or None
    samples = [
        read_sample(args.data, table, bins, "data (linearised)", C_DATA,
                    use_hit_hg=False, mip_cut=args.mip_cut,
                    saturation_adc=args.saturation_adc, max_events=max_data),
        read_sample(args.data, table, bins, "data (raw hit_hg)", C_DATA,
                    use_hit_hg=True, mip_cut=args.mip_cut,
                    saturation_adc=args.saturation_adc, max_events=max_data,
                    dashed=True),
    ]
    for i, (name, path) in enumerate(sim_paths.items()):
        samples.append(
            read_sample(path, table, bins, f"sim ({name})", C_SIM[i % len(C_SIM)],
                        use_hit_hg=False, mip_cut=args.mip_cut,
                        saturation_adc=args.saturation_adc, max_events=None))

    written = [plot_hit_spectrum(samples, bins, args, args.outdir),
               plot_event_sum(samples, bins, args, args.outdir),
               plot_layer_profile(samples, args, args.outdir),
               plot_nhit(samples, bins, args, args.outdir),
               write_layer_csv(samples, args, args.outdir)]

    summary = summarise(samples, bins, args)
    print("\n" + summary)
    summary_path = os.path.join(args.outdir, f"summary_{args.tag}.txt")
    with open(summary_path, "w") as fh:
        fh.write(summary + "\n")
    written.append(summary_path)

    print("\nWritten:")
    for path in written:
        print(f"  {path}")
    return 0


if __name__ == "__main__":
    sys.exit(main())

#!/usr/bin/env python3
"""
Fast shaper against slow shaper, simulation against test beam.

The two shapers do different jobs and are recorded differently on each side:

              fast shaper (trigger)              slow shaper (amplitude)
  simulation  peak [MIP] -> `hit_fast`,          sample at triggerTime+DelayNs
              compared with a fixed Threshold    -> `hit_energy`
  test beam   `hitbit_high`, one BIT per         `adc_high` -> `hit_energy`
              channel and SCA

so the fast shapers cannot be overlaid directly: the test beam keeps only whether
the discriminator fired.  What CAN be compared, and is the thing that matters, is
the **turn-on**: the probability that a cell enters the event at all, as a
function of the amplitude the slow shaper measured for it.

  data  P(hitbit_high = 1 | amplitude)   -- EventBuilder::bestScaPerChannel drops
                                            every channel whose bit never fired,
                                            so this IS the data's hit selection
  sim   P(cell survives the real chain | amplitude), obtained by matching the
        `simple` chain (which keeps every cell, Threshold = 0) against the
        `real` one channel by channel

Everything the simulation does to the fast shaper — the threshold, the noise on
the peak, the 100% efficiency above it — shows up in that curve, and can be put
straight on top of the measured one.

Usage
-----
    python -m analysis.compare_shapers \\
        --sim-real   ecal_sim_e52_real.root \\
        --sim-simple ecal_sim_e52_simple.root \\
        --data-chunk /eos/.../TB2026CERN_run_000013/chunks/chunk_0000.root \\
        --data-events 40 --threshold th230 \\
        --calib-dir ../siwecal-tb2026/calibration/MuonCalib_gaudi \\
        --tag e52 --outdir plots/
"""

from __future__ import annotations

import argparse
import glob
import os
import sys

import numpy as np

N_LAYERS, N_CHIPS, N_CHANS, N_SCAS = 15, 16, 64, 15
MAX_MIP_ADC = 100.0

C_DATA = "#2a78d6"
C_SIM = "#eb6834"
C_SLOW = "#1baf7a"
C_GRID = "#d8d8d4"
C_INK = "#0b0b0b"
C_INK_SOFT = "#52514e"


# --------------------------------------------------------------------------- #
# Calibration
# --------------------------------------------------------------------------- #

def load_tables(calib_dir: str, threshold: str):
    """(MPV[slab,chip,chan], pedestal[slab,chip,chan,sca]) in ADC."""
    mip_files = glob.glob(os.path.join(calib_dir, "mips", threshold, "*highgain*"))
    ped_files = glob.glob(os.path.join(calib_dir, "pedestals", threshold, "*highgain*"))
    if not mip_files or not ped_files:
        sys.exit(f"ERROR: no highgain MIP/pedestal table under {calib_dir} for "
                 f"{threshold}")

    raw = np.loadtxt(mip_files[0], comments="#")
    mpv = np.full((N_LAYERS, N_CHIPS, N_CHANS), np.nan)
    ok = (raw[:, 3] > 0) & (raw[:, 3] <= MAX_MIP_ADC)
    mpv[raw[ok, 0].astype(int), raw[ok, 1].astype(int), raw[ok, 2].astype(int)] = raw[ok, 3]

    rawp = np.loadtxt(ped_files[0], comments="#")
    ped = np.full((N_LAYERS, N_CHIPS, N_CHANS, N_SCAS), np.nan)
    layer, chip, chan = (rawp[:, 0].astype(int), rawp[:, 1].astype(int),
                         rawp[:, 2].astype(int))
    inside = (layer < N_LAYERS) & (chip < N_CHIPS) & (chan < N_CHANS)
    ped[layer[inside], chip[inside], chan[inside], :] = rawp[:, 3::3][inside][:, :N_SCAS]

    print(f"[calib] {threshold}: MPV median {np.nanmedian(mpv):.2f} ADC/MIP, "
          f"pedestal median {np.nanmedian(ped):.1f} ADC")
    return mpv, ped


# --------------------------------------------------------------------------- #
# Data: the measured turn-on
# --------------------------------------------------------------------------- #

def data_turnon(chunks, n_events, mpv, ped, bins_mip, bins_adc=None):
    """Turn-on from the decoded chunks, in MIP and (optionally) in raw ADC.

    One entry per (acquisition, slab, chip, channel): the SCA where that channel
    has the largest pedestal-subtracted ADC, so a signal that also shows up in a
    neighbouring SCA is not counted twice.

    The ADC axis is the calibration-free one — it is the level the discriminator
    actually sits at — while the MIP axis folds in that threshold set's own MPV
    table, which is the axis the reconstruction works on.

    Returns ``{"mip": (num, den), "adc": (num, den) or None, "n": samples}``,
    plus, with ``bins_adc`` given, ``"adc_per_slab": (num, den)`` shaped
    ``(N_LAYERS, nbins)``: the same ADC-axis turn-on slab by slab, for the slab
    that ran at a discriminator of its own (slab 12, the chip-on-board, at
    DAC 243 against the others' 215-230).
    """
    import uproot

    out = {"mip": [np.zeros(len(bins_mip) - 1), np.zeros(len(bins_mip) - 1)],
           "adc": None, "adc_per_slab": None, "n": 0}
    if bins_adc is not None:
        out["adc"] = [np.zeros(len(bins_adc) - 1), np.zeros(len(bins_adc) - 1)]
        out["adc_per_slab"] = [np.zeros((N_LAYERS, len(bins_adc) - 1)),
                               np.zeros((N_LAYERS, len(bins_adc) - 1))]

    ped_t = np.transpose(ped, (0, 1, 3, 2))[None, ...]   # -> (1,slab,chip,sca,chan)
    mpv_b = mpv[None, :, :, None, :]

    for path in chunks:
        arrays = uproot.open(path)["siwecaldecoded"].arrays(
            ["adc_high", "hitbit_high", "badbcid"], entry_stop=n_events,
            library="np")
        adc = np.stack(arrays["adc_high"]).astype(float)
        hitbit = np.stack(arrays["hitbit_high"])
        badbcid = np.stack(arrays["badbcid"])

        pedsub = adc - ped_t
        valid = ((adc > 0) & np.isfinite(pedsub) & (hitbit >= 0)
                 & (badbcid == 0)[..., None])
        pedsub = np.where(valid, pedsub, -np.inf)

        best = np.argmax(pedsub, axis=3)[:, :, :, None, :]
        amp = np.take_along_axis(pedsub, best, axis=3)[:, :, :, 0, :]
        fired = np.take_along_axis(hitbit, best, axis=3)[:, :, :, 0, :]
        seen = np.take_along_axis(valid, best, axis=3)[:, :, :, 0, :]
        mpv_hit = np.take_along_axis(np.broadcast_to(mpv_b, pedsub.shape), best,
                                     axis=3)[:, :, :, 0, :]

        keep = seen & np.isfinite(mpv_hit) & np.isfinite(amp)
        f = fired[keep] == 1
        out["n"] += int(keep.sum())
        slab_of = np.broadcast_to(np.arange(amp.shape[1])[None, :, None, None], amp.shape)[keep]

        if out["adc_per_slab"] is not None:
            nb = len(bins_adc) - 1
            idx = np.digitize(amp[keep], bins_adc) - 1
            inside = (idx >= 0) & (idx < nb)
            flat = slab_of[inside] * nb + idx[inside]
            out["adc_per_slab"][0] += np.bincount(
                flat, weights=f[inside].astype(float), minlength=N_LAYERS * nb).reshape(N_LAYERS, nb)
            out["adc_per_slab"][1] += np.bincount(flat, minlength=N_LAYERS * nb).reshape(N_LAYERS, nb)

        for name, values, edges in (("mip", amp[keep] / mpv_hit[keep], bins_mip),
                                    ("adc", amp[keep], bins_adc)):
            if out[name] is None:
                continue
            idx = np.digitize(values, edges) - 1
            inside = (idx >= 0) & (idx < len(edges) - 1)
            out[name][0] += np.bincount(idx[inside], weights=f[inside].astype(float),
                                        minlength=len(edges) - 1)
            out[name][1] += np.bincount(idx[inside], minlength=len(edges) - 1)

    print(f"[data] {out['n'] / 1e6:.2f}M channel-acquisitions from "
          f"{len(chunks)} chunk(s)")
    return out


# --------------------------------------------------------------------------- #
# Simulation: the turn-on the digitiser applies
# --------------------------------------------------------------------------- #

def _channel_key(slab, chip, chan):
    return (slab.astype(np.int64) * N_CHIPS + chip) * N_CHANS + chan


def sim_turnon(path_real, path_simple, bins, max_events=None):
    """P(cell kept by the real chain | its simple-chain amplitude [MIP]).

    The simple chain runs at Threshold = 0, so it holds every cell with a
    deposit; whatever the real chain is missing was dropped by the fast
    channel's threshold.
    """
    import uproot

    branches = ["nhit_chan", "hit_slab", "hit_chip", "hit_chan", "hit_energy"]
    a_simple = uproot.open(path_simple)["ecal"].arrays(
        branches, entry_stop=max_events, library="np")
    a_real = uproot.open(path_real)["ecal"].arrays(
        branches, entry_stop=max_events, library="np")

    n_events = min(len(a_simple["nhit_chan"]), len(a_real["nhit_chan"]))
    num = np.zeros(len(bins) - 1)
    den = np.zeros(len(bins) - 1)

    for i in range(n_events):
        key_s = _channel_key(a_simple["hit_slab"][i], a_simple["hit_chip"][i],
                             a_simple["hit_chan"][i])
        key_r = _channel_key(a_real["hit_slab"][i], a_real["hit_chip"][i],
                             a_real["hit_chan"][i])
        energy = a_simple["hit_energy"][i].astype(float)
        kept = np.isin(key_s, key_r)

        idx = np.digitize(energy, bins) - 1
        inside = (idx >= 0) & (idx < len(bins) - 1)
        num += np.bincount(idx[inside], weights=kept[inside].astype(float),
                           minlength=len(bins) - 1)
        den += np.bincount(idx[inside], minlength=len(bins) - 1)

    print(f"[sim ] {n_events} events matched channel by channel")
    return num, den


def sim_fast_slow(path_real, max_events=None):
    """The two shapers of the simulation, hit by hit."""
    import uproot

    arrays = uproot.open(path_real)["ecal"].arrays(
        ["hit_energy", "hit_fast", "hit_ismasked"], entry_stop=max_events,
        library="np")
    slow = np.concatenate(arrays["hit_energy"]).astype(float)
    fast = np.concatenate(arrays["hit_fast"]).astype(float)
    masked = np.concatenate(arrays["hit_ismasked"]).astype(int)
    keep = (masked == 0) & (fast > 0)
    if not keep.any():
        sys.exit("ERROR: hit_fast is empty — regenerate the ecal tree with a "
                 "digitisation that writes SiPadHitsRealDigitizedFast")
    return fast[keep], slow[keep]


# --------------------------------------------------------------------------- #
# Plots
# --------------------------------------------------------------------------- #

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


def _efficiency(num, den, min_entries=30):
    eff = np.where(den >= min_entries, num / np.maximum(den, 1), np.nan)
    err = np.where(den >= min_entries,
                   np.sqrt(np.maximum(eff * (1 - eff), 0) / np.maximum(den, 1)),
                   np.nan)
    return eff, err


def plot_turnon(centres, data_eff, data_err, sim_eff, sim_err, args, outdir):
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(8.2, 5.4))
    ax.errorbar(centres, data_eff, yerr=data_err, color=C_DATA, linewidth=2.0,
                marker="o", markersize=3.5, capsize=0,
                label=f"data — hitbit_high ({args.threshold})")
    ax.errorbar(centres, sim_eff, yerr=sim_err, color=C_SIM, linewidth=2.0,
                marker="o", markersize=3.5, capsize=0,
                label=f"sim — fast channel, threshold {args.threshold_mip} MIP")
    ax.axhline(1.0, color=C_INK_SOFT, linewidth=1.0, linestyle=":")
    ax.set_ylim(-0.03, 1.08)
    ax.set_xlim(centres[0], centres[-1])
    _style(ax, "slow-shaper amplitude  [MIP]",
           "P(cell enters the event)",
           f"Trigger turn-on — {args.title}")
    ax.legend(frameon=False, labelcolor=C_INK, loc="lower right")
    fig.tight_layout()
    path = os.path.join(outdir, f"shaper_turnon_{args.tag}.png")
    fig.savefig(path, dpi=140)
    plt.close(fig)
    return path


def plot_fast_slow(fast, slow, args, outdir):
    import matplotlib.pyplot as plt

    fig, (ax, axr, axp) = plt.subplots(1, 3, figsize=(14.0, 4.6))

    hi = float(np.percentile(fast, 99.5))
    ax.hexbin(fast, slow, gridsize=70, extent=(0, hi, 0, hi), bins="log",
              cmap="BuPu", mincnt=1)
    ax.plot([0, hi], [0, hi], color=C_INK_SOFT, linewidth=1.2, linestyle="--")
    _style(ax, "fast-shaper peak  [MIP]",
           "slow-shaper sample  [MIP]", "The two shapers, hit by hit (sim)")
    ax.set_xlim(0, hi)
    ax.set_ylim(0, hi)

    ratio = slow / fast
    axr.hist(ratio, bins=np.linspace(0.8, 1.2, 120), color=C_SIM, histtype="step",
             linewidth=1.8)
    axr.axvline(1.0, color=C_INK_SOFT, linewidth=1.2, linestyle=":")
    _style(axr, "slow / fast", "hits",
           f"median {np.median(ratio):.3f}, IQR "
           f"{np.percentile(ratio, 25):.3f}–{np.percentile(ratio, 75):.3f}")

    # The ratio against amplitude: the hold lands at triggerTime + DelayNs, and a
    # big pulse crosses the threshold sooner, so it is sampled further from the
    # slow peak.  That trend IS the time walk -- the only amplitude dependence
    # the shaping has left now that the kernel is unit-peak normalised.
    edges = np.logspace(np.log10(max(fast.min(), 0.4)), np.log10(hi), 22)
    idx = np.digitize(fast, edges) - 1
    centres, medians = [], []
    for i in range(len(edges) - 1):
        sel = idx == i
        if sel.sum() > 50:
            centres.append(np.sqrt(edges[i] * edges[i + 1]))
            medians.append(float(np.median(ratio[sel])))
    axp.plot(centres, medians, color=C_SIM, linewidth=2.0, marker="o",
             markersize=4)
    axp.axhline(1.0, color=C_INK_SOFT, linewidth=1.2, linestyle=":")
    axp.set_xscale("log")
    _style(axp, "fast-shaper peak  [MIP]", "median slow / fast",
           f"Time walk: hold at {args.delay_ns:.0f} ns, slow peak at "
           f"{args.tau_slow_ns:.0f} ns")

    fig.tight_layout()
    path = os.path.join(outdir, f"shaper_fast_vs_slow_{args.tag}.png")
    fig.savefig(path, dpi=140)
    plt.close(fig)
    return path


def plot_spectra(fast, slow, args, outdir):
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(8.2, 5.2))
    bins = np.linspace(0, float(np.percentile(slow, 99.5)), 160)
    ax.hist(fast, bins=bins, histtype="step", linewidth=1.8, color=C_SIM,
            label="fast channel (peak)")
    ax.hist(slow, bins=bins, histtype="step", linewidth=1.8, color=C_SLOW,
            label="slow channel (sample)")
    ax.axvline(args.threshold_mip, color=C_INK_SOFT, linewidth=1.2, linestyle=":")
    ax.set_yscale("log")
    ax.set_xlim(bins[0], bins[-1])
    _style(ax, "amplitude  [MIP]", "hits",
           f"Both channels' spectra (sim) — {args.title}")
    ax.legend(frameon=False, labelcolor=C_INK)
    fig.tight_layout()
    path = os.path.join(outdir, f"shaper_spectra_{args.tag}.png")
    fig.savefig(path, dpi=140)
    plt.close(fig)
    return path


# --------------------------------------------------------------------------- #
# CLI
# --------------------------------------------------------------------------- #

def _crossing(centres, eff, level):
    """First bin where the efficiency reaches `level`, interpolated."""
    finite = np.isfinite(eff)
    idx = np.where(finite & (eff >= level))[0]
    if not len(idx):
        return float("nan")
    i = idx[0]
    if i == 0 or not np.isfinite(eff[i - 1]) or eff[i] == eff[i - 1]:
        return float(centres[i])
    frac = (level - eff[i - 1]) / (eff[i] - eff[i - 1])
    return float(centres[i - 1] + frac * (centres[i] - centres[i - 1]))


def main(argv=None) -> int:
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--sim-real", required=True, help="ecal tree of the real chain")
    p.add_argument("--sim-simple", required=True,
                   help="ecal tree of the simple chain (same sample, Threshold=0)")
    p.add_argument("--data-chunk", action="append", default=[], required=True,
                   help="Decoded chunk (siwecaldecoded tree), repeatable")
    p.add_argument("--data-events", type=int, default=40,
                   help="Acquisitions read per chunk (default: 40)")
    p.add_argument("--calib-dir",
                   default="../siwecal-tb2026/calibration/MuonCalib_gaudi")
    p.add_argument("--threshold", default="th230",
                   help="Calibration set matching the data run (default: th230)")
    p.add_argument("--threshold-mip", type=float, default=0.5,
                   help="The simulation's fast-channel threshold [MIP]")
    p.add_argument("--max-sim-events", type=int, default=None)
    p.add_argument("--delay-ns", type=float, default=160.0,
                   help="RealDigitizer's DelayNs, for the time-walk panel")
    p.add_argument("--tau-slow-ns", type=float, default=180.0,
                   help="RealDigitizer's TauSlowNs, for the time-walk panel")
    p.add_argument("--tag", default="shapers")
    p.add_argument("--title", default="")
    p.add_argument("--outdir", default="plots")
    args = p.parse_args(argv)
    if not args.title:
        args.title = args.tag

    import matplotlib
    matplotlib.use("Agg")

    os.makedirs(args.outdir, exist_ok=True)
    mpv, ped = load_tables(args.calib_dir, args.threshold)

    bins = np.arange(0.0, 5.001, 0.1)
    centres = 0.5 * (bins[:-1] + bins[1:])

    d_num, d_den = data_turnon(args.data_chunk, args.data_events, mpv, ped,
                               bins)["mip"]
    s_num, s_den = sim_turnon(args.sim_real, args.sim_simple, bins,
                              args.max_sim_events)
    data_eff, data_err = _efficiency(d_num, d_den)
    sim_eff, sim_err = _efficiency(s_num, s_den)

    fast, slow = sim_fast_slow(args.sim_real, args.max_sim_events)

    written = [plot_turnon(centres, data_eff, data_err, sim_eff, sim_err,
                           args, args.outdir),
               plot_fast_slow(fast, slow, args, args.outdir),
               plot_spectra(fast, slow, args, args.outdir)]

    plateau = np.isfinite(data_eff) & (centres > 2.5)
    sim_plateau = np.isfinite(sim_eff) & (centres > 2.5)
    lines = [
        f"# Fast vs slow shaper — {args.title}",
        f"#   data: {len(args.data_chunk)} chunk(s) x {args.data_events} "
        f"acquisitions, {args.threshold}",
        f"#   sim : {os.path.basename(args.sim_real)} vs "
        f"{os.path.basename(args.sim_simple)}, threshold {args.threshold_mip} MIP",
        "",
        f"{'quantity':<34}{'data':>10}{'sim':>10}",
        f"{'turn-on 50% [MIP]':<34}{_crossing(centres, data_eff, 0.5):>10.2f}"
        f"{_crossing(centres, sim_eff, 0.5):>10.2f}",
        f"{'turn-on 90% [MIP]':<34}{_crossing(centres, data_eff, 0.9):>10.2f}"
        f"{_crossing(centres, sim_eff, 0.9):>10.2f}",
        f"{'plateau efficiency (>2.5 MIP)':<34}"
        f"{np.nanmean(data_eff[plateau]):>10.3f}"
        f"{np.nanmean(sim_eff[sim_plateau]):>10.3f}",
        "",
        f"{'sim slow/fast ratio, median':<34}{'':>10}{np.median(slow / fast):>10.3f}",
        f"{'sim slow/fast ratio, IQR':<34}{'':>10}"
        f"{np.percentile(slow / fast, 75) - np.percentile(slow / fast, 25):>10.3f}",
    ]
    summary = "\n".join(lines)
    print("\n" + summary)
    summary_path = os.path.join(args.outdir, f"summary_shapers_{args.tag}.txt")
    with open(summary_path, "w") as fh:
        fh.write(summary + "\n")
    written.append(summary_path)

    print("\nWritten:")
    for path in written:
        print(f"  {path}")
    return 0


if __name__ == "__main__":
    sys.exit(main())

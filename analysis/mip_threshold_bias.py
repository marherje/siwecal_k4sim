#!/usr/bin/env python3
"""
How much does the trigger threshold bias the MIP calibration?

A MIP calibration fits (or peak-finds) the MPV of the muon ADC spectrum.  When
the trigger threshold sits near or above that peak, the spectrum the calibration
sees is the Landau with its rising edge cut off, and the recovered MPV runs away
with the threshold instead of measuring the peak.  Every hit energy in that run
is then divided by an MPV that is too large, so the whole energy scale of the
data comes out low — flat in beam energy, which is exactly the signature that a
saturation effect does NOT have.

This script measures that bias, using the *simulation* as the ground truth: a
simulated MIP-like sample has a known MIP scale (1 MIP by construction, see
analysis/tests/test_cell_shaping.py), so pushing it through a threshold and
recovering the MPV the way the calibration does gives the bias directly.

The estimator here is PedestalMipCalibrator's fallback: the maximum of the
spectrum rebinned to ``--rebin`` ADC (status 2, ``empv = -2`` in the table), which
is what a large fraction of the high-threshold channels end up using when the
per-channel Landau-Gauss fit is not believed.

Usage
-----
    python -m analysis.mip_threshold_bias \\
        --sim  ecal_sim_mu100_simple.root \\
        --mip-file ../siwecal-tb2026/calibration/MuonCalib_gaudi/mips/th220/\\
MIP_pedestalsubmode1_TB2026CERN_run_000th220_highgain.txt \\
        --target-bias 1.47 --outdir plots/
"""

from __future__ import annotations

import argparse
import os
import sys

import numpy as np

N_LAYERS, N_CHIPS, N_CHANS = 15, 16, 64
MAX_MIP_ADC = 100.0

C_LINE = "#eb6834"
C_REF = "#2a78d6"
C_GRID = "#d8d8d4"
C_INK = "#0b0b0b"
C_INK_SOFT = "#52514e"


def load_mip_table(path: str) -> np.ndarray:
    raw = np.loadtxt(path, comments="#")
    table = np.full((N_LAYERS, N_CHIPS, N_CHANS), np.nan)
    layer, chip, chan, mpv = (raw[:, 0].astype(int), raw[:, 1].astype(int),
                              raw[:, 2].astype(int), raw[:, 3].astype(float))
    ok = (mpv > 0) & (mpv <= MAX_MIP_ADC) & np.isfinite(mpv)
    ok &= (layer >= 0) & (layer < N_LAYERS) & (chip >= 0) & (chip < N_CHIPS)
    ok &= (chan >= 0) & (chan < N_CHANS)
    table[layer[ok], chip[ok], chan[ok]] = mpv[ok]
    return table


def peak_estimate(values: np.ndarray, rebin: float, hi: float) -> float:
    """The calibration's fallback MPV: the maximum of the rebinned spectrum."""
    bins = np.arange(0.0, hi + rebin, rebin)
    counts, _ = np.histogram(values, bins=bins)
    if counts.sum() == 0:
        return float("nan")
    i = int(np.argmax(counts))
    return 0.5 * (bins[i] + bins[i + 1])


def main(argv=None) -> int:
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--sim", required=True,
                   help="Simulated MIP-like ecal tree (e.g. the muon sample)")
    p.add_argument("--mip-file", required=True,
                   help="MIP table that defines the ADC per MIP of each channel")
    p.add_argument("--pedestal-noise", type=float, default=2.0,
                   help="Gaussian pedestal noise added in ADC (default: 2.0, "
                        "the measured median pedestal width)")
    p.add_argument("--rebin", type=float, default=4.0,
                   help="ADC per bin of the peak estimator (PedestalMipCalibrator's "
                        "MipPeakRebin, default 4)")
    p.add_argument("--target-bias", type=float, default=None,
                   help="Report the threshold that would produce this MPV bias "
                        "(e.g. 1.47, the measured th230/th220 ratio)")
    p.add_argument("--seed", type=int, default=7)
    p.add_argument("--outdir", default="")
    p.add_argument("--tag", default="mip_threshold_bias")
    args = p.parse_args(argv)

    import uproot

    table = load_mip_table(args.mip_file)
    arrays = uproot.open(args.sim)["ecal"].arrays(
        ["hit_slab", "hit_chip", "hit_chan", "hit_energy", "hit_ismasked"],
        library="np")
    slab = np.concatenate(arrays["hit_slab"]).astype(int)
    chip = np.concatenate(arrays["hit_chip"]).astype(int)
    chan = np.concatenate(arrays["hit_chan"]).astype(int)
    energy = np.concatenate(arrays["hit_energy"]).astype(float)
    masked = np.concatenate(arrays["hit_ismasked"]).astype(int)

    ok = ((slab >= 0) & (slab < N_LAYERS) & (chip >= 0) & (chip < N_CHIPS)
          & (chan >= 0) & (chan < N_CHANS) & (masked == 0))
    mpv = np.full(len(slab), np.nan)
    mpv[ok] = table[slab[ok], chip[ok], chan[ok]]
    keep = ok & np.isfinite(mpv) & (energy > 0)

    rng = np.random.default_rng(args.seed)
    adc = energy[keep] * mpv[keep] + rng.normal(0.0, args.pedestal_noise, int(keep.sum()))
    true_mpv = float(np.nanmedian(mpv[keep]))
    hi = 6.0 * true_mpv

    thresholds = np.concatenate((np.arange(0.0, 1.0, 0.1),
                                 np.arange(1.0, 2.01, 0.05))) * true_mpv
    rows = []
    for thr in thresholds:
        sel = adc > thr
        if sel.sum() < 100:
            continue
        rows.append((thr, thr / true_mpv, sel.mean(),
                     peak_estimate(adc[sel], args.rebin, hi)))

    print(f"[bias] {len(adc)} simulated MIP hits, true MPV = {true_mpv:.2f} ADC "
          f"(median of {os.path.basename(args.mip_file)})")
    print(f"{'thr[ADC]':>9}{'thr/MPV':>9}{'kept':>8}{'peak[ADC]':>11}{'bias':>8}")
    for thr, ratio, kept, peak in rows:
        print(f"{thr:>9.1f}{ratio:>9.2f}{kept * 100:>7.1f}%{peak:>11.2f}"
              f"{peak / true_mpv:>8.2f}")

    if args.target_bias:
        bias = np.array([r[3] / true_mpv for r in rows])
        ratio = np.array([r[1] for r in rows])
        reached = np.where(bias >= args.target_bias)[0]
        if len(reached):
            print(f"\n[bias] a bias of x{args.target_bias:.2f} needs a threshold of "
                  f"{ratio[reached[0]]:.2f} x MPV "
                  f"({ratio[reached[0]] * true_mpv:.1f} ADC), which keeps "
                  f"{rows[reached[0]][2] * 100:.0f}% of the MIP hits")
        else:
            print(f"\n[bias] x{args.target_bias:.2f} is not reached below "
                  f"{ratio[-1]:.2f} x MPV")

    if args.outdir:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        os.makedirs(args.outdir, exist_ok=True)
        fig, ax = plt.subplots(figsize=(8.2, 5.0))
        ax.plot([r[1] for r in rows], [r[3] / true_mpv for r in rows],
                color=C_LINE, linewidth=2.0, label="recovered MPV / true MPV")
        ax.axhline(1.0, color=C_INK_SOFT, linewidth=1.0, linestyle=":")
        if args.target_bias:
            ax.axhline(args.target_bias, color=C_REF, linewidth=1.6,
                       linestyle="--",
                       label=f"measured th230/th220 = {args.target_bias:.2f}")
        ax.set_xlabel("trigger threshold  [units of the true MIP MPV]", color=C_INK)
        ax.set_ylabel("MPV recovered / true MPV", color=C_INK)
        ax.set_title("A threshold above the MIP peak inflates the calibration MPV",
                     color=C_INK, fontsize=11, pad=10)
        ax.grid(True, color=C_GRID, linewidth=0.6)
        ax.set_axisbelow(True)
        for spine in ("top", "right"):
            ax.spines[spine].set_visible(False)
        for spine in ("left", "bottom"):
            ax.spines[spine].set_color(C_GRID)
        ax.tick_params(colors=C_INK_SOFT)
        ax.legend(frameon=False, labelcolor=C_INK)
        fig.tight_layout()
        out = os.path.join(args.outdir, f"{args.tag}.png")
        fig.savefig(out, dpi=140)
        plt.close(fig)
        print(f"\nWritten: {out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())

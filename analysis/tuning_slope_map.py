#!/usr/bin/env python
"""
sim_tuning 2.1: per-chip baseline slope map dB/dQ for the emulation (sim_tuning.apply_baseline slope_map), from the
records of analysis/cm_trigger_test.py (retrigger columns vetoed), the available runs pooled.

    python analysis/tuning_slope_map.py <cm npz dir> <out.npz> [runs, default e291,e287,e286,e166]

Writes slope[slab, chip] in ADC per ADC (NaN where the chip never reaches Q > 1500 ADC often enough), its error, the
pooled all-chip slope, and slope_map.png next to the npz.
"""
import os
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from analysis.cm_per_chip import load, slope  # noqa: E402


def main():
    indir, out = sys.argv[1], sys.argv[2]
    runs = (sys.argv[3] if len(sys.argv) > 3 else "e291,e287,e286,e166").split(",")
    runs = [r for r in runs if os.path.exists(f"{indir}/cm_trigger_{r}.npz")]
    D = [load(indir, r) for r in runs]
    P = {k: np.concatenate([d[k] for d in D]) for k in D[0]}
    S = np.full((15, 16), np.nan); E = np.full((15, 16), np.nan)
    for s in range(15):
        for c in range(16):
            m = (P["slab"] == s) & (P["chip"] == c)
            if m.sum() >= 200:
                S[s, c], E[s, c] = slope(P["q"][m], P["base_far"][m])
    g, ge = slope(P["q"], P["base_far"], nmin=200)
    np.savez(out, slope=S / 1000.0, err=E / 1000.0, pooled=g / 1000.0, runs=np.array(runs))
    print(f"runs {runs}: pooled {g:.3f} +- {ge:.3f} ADC per 1000 ADC; {np.isfinite(S).sum()} chips measured, "
          f"median {np.nanmedian(S):.3f}, 16-84 % {np.nanpercentile(S, 16):.3f} .. {np.nanpercentile(S, 84):.3f}")
    fig, ax = plt.subplots(figsize=(7, 4.5))
    im = ax.imshow(S.T, origin="lower", cmap="viridis_r", vmin=-3, vmax=0, aspect="auto")
    fig.colorbar(im, ax=ax, label="dB/dQ [ADC per 1000 ADC]")
    ax.set_xlabel("slab"); ax.set_ylabel("chip"); ax.set_title(f"chip baseline slope ({', '.join(runs)})")
    fig.tight_layout(); fig.savefig(os.path.join(os.path.dirname(out), "slope_map.png"), dpi=120); plt.close(fig)


if __name__ == "__main__":
    main()

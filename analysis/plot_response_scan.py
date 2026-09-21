#!/usr/bin/env python3
"""
Calorimeter response (MIP per GeV) versus beam energy, data against simulation.

This is the plot that separates a *scale* problem from a *saturation* problem:
a flat line that sits too low is a calibration / sampling-response question,
while a line that bends over at high energy is the detector or the event
building losing signal where the occupancy is highest.

Both sides are read in MIP (``sum_energy``), so nothing here depends on the
MIP -> ADC conversion that analysis/compare_adc_data_sim.py needs; the only
assumption is that the simulation's MIP and the data's MIP are the same unit,
which the muon comparison in that script is there to check.

Usage
-----
    python -m analysis.plot_response_scan \\
        --data 7.5=/eos/.../ecal_TB2026CERN_run_000022.root \\
        --data 52=/eos/.../ecal_TB2026CERN_run_000013.root \\
        --sim  52:digi=/path/ecal_sim_e52_real.root \\
        --outdir plots/ --tag P1
"""

from __future__ import annotations

import argparse
import os
import sys
from typing import Dict, List, Tuple

import numpy as np

C_DATA = "#2a78d6"
C_SIM = ["#eb6834", "#1baf7a", "#eda100"]
C_GRID = "#d8d8d4"
C_INK = "#0b0b0b"
C_INK_SOFT = "#52514e"


def _read_mean(path: str, max_events: int) -> Tuple[float, float, int]:
    """Mean and RMS of the per-event sum_energy [MIP], and the event count."""
    import uproot

    tree = uproot.open(path)["ecal"]
    stop = min(tree.num_entries, max_events) if max_events else tree.num_entries
    values = tree["sum_energy"].array(entry_stop=stop, library="np")
    return float(values.mean()), float(values.std()), len(values)


def _parse_points(specs: List[str], with_chain: bool) -> List[Tuple]:
    out = []
    for spec in specs:
        if "=" not in spec:
            sys.exit(f"ERROR: expected ENERGY{'[:CHAIN]' if with_chain else ''}"
                     f"=PATH, got '{spec}'")
        key, path = spec.split("=", 1)
        chain = ""
        if with_chain and ":" in key:
            key, chain = key.split(":", 1)
        out.append((float(key), chain, path))
    return sorted(out)


def main(argv=None) -> int:
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--data", action="append", default=[], metavar="E=PATH",
                   required=True, help="Data ecal tree at beam energy E [GeV]")
    p.add_argument("--sim", action="append", default=[], metavar="E:CHAIN=PATH",
                   help="Simulation ecal tree; CHAIN names the legend entry")
    p.add_argument("--max-events", type=int, default=30000,
                   help="Events read per file (0 = all; default: 30000)")
    p.add_argument("--tag", default="scan")
    p.add_argument("--title", default="")
    p.add_argument("--outdir", default="plots")
    args = p.parse_args(argv)

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    os.makedirs(args.outdir, exist_ok=True)

    data_points = _parse_points(args.data, with_chain=False)
    sim_points = _parse_points(args.sim, with_chain=True)

    rows = []
    series: Dict[str, List[Tuple[float, float, float, int]]] = {}
    for energy, _, path in data_points:
        mean, rms, n = _read_mean(path, args.max_events)
        series.setdefault("data", []).append((energy, mean, rms, n))
        rows.append(("data", energy, mean, rms, n, os.path.basename(path)))
    for energy, chain, path in sim_points:
        mean, rms, n = _read_mean(path, args.max_events)
        name = f"sim ({chain})" if chain else "sim"
        series.setdefault(name, []).append((energy, mean, rms, n))
        rows.append((name, energy, mean, rms, n, os.path.basename(path)))

    fig, (ax, axr) = plt.subplots(2, 1, figsize=(8.2, 7.0), sharex=True,
                                  gridspec_kw={"height_ratios": [2, 1]})
    colors = {"data": C_DATA}
    for i, name in enumerate(n for n in series if n != "data"):
        colors[name] = C_SIM[i % len(C_SIM)]

    for name, points in series.items():
        e = np.array([q[0] for q in points])
        mu = np.array([q[1] for q in points])
        rms = np.array([q[2] for q in points])
        nev = np.array([q[3] for q in points], dtype=float)
        order = np.argsort(e)
        e, mu, rms, nev = e[order], mu[order], rms[order], nev[order]
        # Error on the MEAN response: the spread of the sample is the bottom
        # panel's job, and drawing it here would read as an uncertainty.
        ax.errorbar(e, mu / e, yerr=rms / e / np.sqrt(np.maximum(nev, 1)),
                    color=colors[name], linewidth=2.0, marker="o", markersize=5,
                    capsize=0, label=name)
        axr.plot(e, rms / mu, color=colors[name], linewidth=2.0, marker="o",
                 markersize=5)

    for axis, xlabel, ylabel in ((ax, "", "response  [MIP / GeV]"),
                                 (axr, "beam energy  [GeV]", "sigma / mu")):
        axis.set_xlabel(xlabel, color=C_INK)
        axis.set_ylabel(ylabel, color=C_INK)
        axis.grid(True, color=C_GRID, linewidth=0.6)
        axis.set_axisbelow(True)
        for spine in ("top", "right"):
            axis.spines[spine].set_visible(False)
        for spine in ("left", "bottom"):
            axis.spines[spine].set_color(C_GRID)
        axis.tick_params(colors=C_INK_SOFT)
    ax.set_ylim(bottom=0)
    ax.set_title(args.title or "Response vs beam energy", color=C_INK,
                 fontsize=11, pad=10)
    ax.legend(frameon=False, labelcolor=C_INK)
    fig.tight_layout()
    path = os.path.join(args.outdir, f"response_scan_{args.tag}.png")
    fig.savefig(path, dpi=140)
    plt.close(fig)

    lines = [f"{'sample':<14}{'E[GeV]':>8}{'<sum>[MIP]':>12}{'MIP/GeV':>10}"
             f"{'sigma/mu':>10}{'events':>9}  file"]
    for name, energy, mean, rms, n, base in rows:
        lines.append(f"{name:<14}{energy:>8.1f}{mean:>12.1f}{mean / energy:>10.1f}"
                     f"{rms / mean:>10.3f}{n:>9d}  {base}")
    summary = "\n".join(lines)
    print(summary)
    with open(os.path.join(args.outdir, f"response_scan_{args.tag}.txt"), "w") as fh:
        fh.write(summary + "\n")
    print(f"\nWritten: {path}")
    return 0


if __name__ == "__main__":
    sys.exit(main())

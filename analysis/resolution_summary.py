#!/usr/bin/env python3
"""
Resolution and response against beam energy, data and simulation, from the
`summary_event_<tag>.txt` files written by `analysis/compare_event_level.py`.

    python -m analysis.resolution_summary --outdir plots/ \\
        --point th230_e20=20 --point th230_e52=52 --point th220_e74=74 ...

Writes resolution_vs_energy_<tag>.png (sigma/mu of the Gaussian core for the
ADC sum, the MIP sum and the hit count) and response_vs_energy_<tag>.png
(core mu per GeV), plus a markdown table.
"""

from __future__ import annotations

import argparse
import os
import re
import sys

import numpy as np

C_DATA = "#2a78d6"
C_SIM = "#d2551f"
MARK = {"th210": "o", "th220": "s", "th230": "^"}


def parse(path):
    rows = {}
    for line in open(path):
        m = re.match(r"^(ADC|MIP|hits)\s+(data|sim)\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)", line)
        if m:
            q, side = m.group(1), m.group(2)
            rows[(q, side)] = dict(mean=float(m.group(3)), peak=float(m.group(4)),
                                   mu=float(m.group(5)), sigma=float(m.group(6)),
                                   res=float(m.group(7)), raw=float(m.group(8)))
    return rows


def main(argv=None) -> int:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--point", action="append", default=[], metavar="TAG=E_GeV", required=True)
    p.add_argument("--outdir", default="plots")
    p.add_argument("--tag", default="final")
    args = p.parse_args(argv)

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    points = []
    for spec in args.point:
        tag, _, e = spec.partition("=")
        rows = parse(os.path.join(args.outdir, f"summary_event_{tag}.txt"))
        th = tag.split("_")[0]
        points.append((tag, th, float(e), rows))

    fig, axes = plt.subplots(1, 3, figsize=(15, 4.6))
    for ax, q, label in zip(axes, ("ADC", "MIP", "hits"),
                            ("event ADC sum", "event energy sum [MIP]", "hits per event")):
        for side, colour in (("data", C_DATA), ("sim", C_SIM)):
            for tag, th, e, rows in points:
                r = rows[(q, side)]
                ax.plot(e, r["res"], MARK[th], color=colour, markersize=8,
                        label=f"{side} {th}" if True else None)
        ax.set_xlabel("beam energy [GeV]")
        ax.set_ylabel("σ/μ of the Gaussian core")
        ax.set_title(f"resolution — {label}")
        ax.set_ylim(0, None)
        ax.grid(alpha=0.3)
        h, l = ax.get_legend_handles_labels()
        seen = {}
        for hh, ll in zip(h, l):
            seen.setdefault(ll, hh)
        ax.legend(seen.values(), seen.keys(), fontsize=7.5, ncol=2)
    fig.tight_layout()
    out1 = os.path.join(args.outdir, f"resolution_vs_energy_{args.tag}.png")
    fig.savefig(out1, dpi=120)
    plt.close(fig)

    fig, axes = plt.subplots(1, 2, figsize=(11, 4.6))
    for ax, q, label in zip(axes, ("MIP", "ADC"), ("MIP / GeV", "ADC / GeV")):
        for side, colour in (("data", C_DATA), ("sim", C_SIM)):
            for tag, th, e, rows in points:
                ax.plot(e, rows[(q, side)]["mu"] / e, MARK[th], color=colour, markersize=8,
                        label=f"{side} {th}")
        ax.set_xlabel("beam energy [GeV]")
        ax.set_ylabel(f"response, core μ per GeV [{label}]")
        ax.set_ylim(0, None)
        ax.grid(alpha=0.3)
        h, l = ax.get_legend_handles_labels()
        seen = {}
        for hh, ll in zip(h, l):
            seen.setdefault(ll, hh)
        ax.legend(seen.values(), seen.keys(), fontsize=7.5, ncol=2)
    fig.tight_layout()
    out2 = os.path.join(args.outdir, f"response_vs_energy_{args.tag}.png")
    fig.savefig(out2, dpi=120)
    plt.close(fig)

    lines = ["| point | E [GeV] | ADC μ d / s | sim/d | σ/μ ADC d / s | MIP μ d / s | sim/d | σ/μ MIP d / s | hits d / s | sim/d | σ/μ hits d / s |",
             "|---|---|---|---|---|---|---|---|---|---|---|"]
    for tag, th, e, r in sorted(points, key=lambda t: (t[1], t[2])):
        f = lambda q, k: (r[(q, "data")][k], r[(q, "sim")][k])  # noqa: E731
        lines.append(f"| {tag} | {e:g} | {f('ADC','mu')[0]:.0f} / {f('ADC','mu')[1]:.0f} | {f('ADC','mu')[1]/f('ADC','mu')[0]:.3f} "
                     f"| {f('ADC','res')[0]:.3f} / {f('ADC','res')[1]:.3f} "
                     f"| {f('MIP','mu')[0]:.0f} / {f('MIP','mu')[1]:.0f} | {f('MIP','mu')[1]/f('MIP','mu')[0]:.3f} "
                     f"| {f('MIP','res')[0]:.3f} / {f('MIP','res')[1]:.3f} "
                     f"| {f('hits','mu')[0]:.0f} / {f('hits','mu')[1]:.0f} | {f('hits','mu')[1]/f('hits','mu')[0]:.3f} "
                     f"| {f('hits','res')[0]:.3f} / {f('hits','res')[1]:.3f} |")
    text = "\n".join(lines)
    print(text)
    with open(os.path.join(args.outdir, f"resolution_table_{args.tag}.md"), "w") as fh:
        fh.write(text + "\n")
    print(f"Written: {out1}\n         {out2}")
    return 0


if __name__ == "__main__":
    sys.exit(main())

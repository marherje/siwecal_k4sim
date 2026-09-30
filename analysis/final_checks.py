#!/usr/bin/env python
"""
Closing checks of the final sample (from the valtrees):

  positions   per group: median shower barycentre (bar_x, bar_y) of the data (Reconstructed_final, pooled runs) and of
              the digitised simulation, against the nominal beam centre of the sim sample -> positions.txt
  hitsel      per run: hit-bit (Reconstructed_final) against adc (Reconstructed_final_adc) selection, shower-like
              events: Gaussian-core energy [MIP], sigma/mu, median nhit -> hitsel.txt, hitsel.png

    python analysis/final_checks.py <outdir> [positions] [hitsel]
"""
import os
import sys

import numpy as np
import uproot
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(_HERE))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(_HERE)), "siwecal-tb2026"))
from analysis.compare_event_level import gauss_core  # noqa: E402
from analysis.final_panels import load_groups, sim_valtree  # noqa: E402

B = "/eos/experiment/drdcalo/siw-ecal/TB2026-06"
RECO, RECO_ADC = f"{B}/Reconstructed_final", f"{B}/Reconstructed_final_adc"


def read(paths, branches):
    cols = {b: [] for b in branches}
    for p in paths:
        for a in uproot.iterate(f"{p}:ecal", branches, step_size="400 MB", library="np"):
            for b in branches:
                cols[b].append(a[b])
    c = {b: np.concatenate(v).astype(float) for b, v in cols.items()}
    k = c["nhit_chan"] > 0.5 * np.percentile(c["nhit_chan"], 90)
    return {b: v[k] for b, v in c.items()}


def positions(out):
    lines = [f"{'group':<24}{'sim beam x,y':>14}{'data bar x,y':>16}{'sim bar x,y':>16}{'|data-sim| [mm]':>17}"]
    worst = 0.0
    for g in load_groups():
        vt = [f"{RECO}/{r}/ecal_{r}.valtree.root" for r in g["runs"]]
        d = read(vt, ["nhit_chan", "bar_x", "bar_y"])
        s = read([sim_valtree(g, "digi")], ["nhit_chan", "bar_x", "bar_y"])
        dx, dy = np.median(d["bar_x"]), np.median(d["bar_y"])
        sx, sy = np.median(s["bar_x"]), np.median(s["bar_y"])
        dist = float(np.hypot(dx - sx, dy - sy)); worst = max(worst, dist)
        bx, by = g["pos"].split("_")
        lines.append(f"{g['name']:<24}{bx + ',' + by:>14}{dx:8.1f},{dy:6.1f} {sx:8.1f},{sy:6.1f}{dist:15.1f}")
        print(lines[-1], flush=True)
    lines.append(f"\nlargest data-sim barycentre distance: {worst:.1f} mm")
    open(f"{out}/positions.txt", "w").write("\n".join(lines) + "\n")
    print(lines[-1])


def hitsel(out):
    runs = sorted({r for g in load_groups() for r in g["runs"]})
    E = {}
    for g in load_groups():
        for r in g["runs"]:
            E[r] = float(g["E"])
    rows, lines = [], [f"{'run':<24}{'E':>5}{'mu hitbit':>11}{'mu adc':>9}{'adc/hb':>8}{'res hb':>8}{'res adc':>8}"
                            f"{'nhit hb':>9}{'nhit adc':>9}{'adc/hb':>8}"]
    for r in runs:
        h = read([f"{RECO}/{r}/ecal_{r}.valtree.root"], ["nhit_chan", "sum_energy"])
        a = read([f"{RECO_ADC}/{r}/ecal_{r}.valtree.root"], ["nhit_chan", "sum_energy"])
        mh, sh = gauss_core(h["sum_energy"]); ma, sa = gauss_core(a["sum_energy"])
        nh, na = np.median(h["nhit_chan"]), np.median(a["nhit_chan"])
        rows.append((E[r], ma / mh, na / nh, sh / mh, sa / ma, r))
        lines.append(f"{r.replace('TB2026CERN_', ''):<24}{E[r]:5g}{mh:11.0f}{ma:9.0f}{ma / mh:8.3f}{sh / mh:8.3f}"
                     f"{sa / ma:8.3f}{nh:9.0f}{na:9.0f}{na / nh:8.3f}")
        print(lines[-1], flush=True)
    rr = np.array([x[1] for x in rows]); nn = np.array([x[2] for x in rows])
    lines.append(f"\nadc / hit-bit energy: median {np.median(rr):.3f}, range {rr.min():.3f}-{rr.max():.3f}; "
                 f"nhit: median {np.median(nn):.3f}, range {nn.min():.3f}-{nn.max():.3f}")
    open(f"{out}/hitsel.txt", "w").write("\n".join(lines) + "\n")
    print(lines[-1])
    fig, ax = plt.subplots(1, 2, figsize=(13.5, 5))
    e = np.array([x[0] for x in rows])
    ax[0].plot(e, rr, "o", color="#29a3dc", ms=8, label="energy [MIP], Gaussian core")
    ax[0].plot(e * 1.03, nn, "s", color="#c0392b", ms=7, mfc="none", mew=1.8, label="number of hits, median")
    ax[0].axhspan(1.04, 1.08, color="0.5", alpha=0.15, label="expected from the retrigger (4-8 %)")
    ax[0].axhline(1, color="0.5", ls="--", lw=1)
    ax[0].set_ylabel("adc selection / hit-bit selection"); ax[0].legend(fontsize=9)
    ax[1].plot(e, [100 * x[3] for x in rows], "o", color="#29a3dc", ms=8, label="hit-bit selection")
    ax[1].plot(e * 1.03, [100 * x[4] for x in rows], "s", color="#c0392b", ms=7, mfc="none", mew=1.8, label="adc selection")
    ax[1].set_ylabel("energy resolution σ/μ [%]"); ax[1].legend(fontsize=9)
    for a_ in ax:
        a_.set_xscale("log"); a_.set_xticks([7.5, 10, 20, 34, 52, 74, 99]); a_.set_xticklabels(["7.5", "10", "20", "34", "52", "74", "99"])
        a_.set_xlabel("Beam energy [GeV]"); a_.grid(alpha=0.25, which="both")
    fig.suptitle("Reconstructed_final (hit bit) vs Reconstructed_final_adc (hit bit or ADC − pedestal > 30), every run, shower-like events",
                 fontsize=11.5)
    fig.tight_layout(rect=[0, 0, 1, 0.93]); fig.savefig(f"{out}/hitsel.png", dpi=110)


if __name__ == "__main__":
    out = sys.argv[1]; os.makedirs(out, exist_ok=True)
    todo = sys.argv[2:] or ["positions", "hitsel"]
    if "positions" in todo:
        positions(out)
    if "hitsel" in todo:
        hitsel(out)

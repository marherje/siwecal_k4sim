#!/usr/bin/env python
"""
Is there a class of data events that looks like the simulation? Shower variables of data events split by how they
were recorded (events in the same DAQ acquisition, order in it, time in the run) and by where they hit, each against
the digitised simulation (final_v6). Position classes cut the simulation the same way; the other classes are
compared to the whole simulation.

    python analysis/event_classes.py <outdir> --group th210_S1_e20 --group th210_S1_e74

Data: Reconstruction ecal tree (n_events_acq, spill, bcid) + valtree (shower variables, 0.5 MIP cut), same
entries. Selection: final_panels.py's shower-like cut. Writes event_classes.png and event_classes.txt.
"""
import argparse
import os
import sys

import numpy as np
import uproot
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(_HERE))
os.environ.setdefault("FINAL_GROUPS", os.path.join(_HERE, "final_v6_groups.txt"))
os.environ.setdefault("FINAL_SIM", "/eos/experiment/drdcalo/siw-ecal/TB2026-06/Simulation/Processed/adc_vs_tb/final_v6")
from analysis.final_panels import RECO, load_groups, sim_valtree, showerlike  # noqa: E402
from analysis.compare_event_level import gauss_core  # noqa: E402

INK, MUTED, GRID = "#1b2227", "#56626b", "#d7dde2"
plt.rcParams.update({"font.size": 11, "axes.edgecolor": MUTED, "axes.labelcolor": INK, "xtick.color": MUTED,
                     "ytick.color": MUTED, "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.7,
                     "axes.spines.top": False, "axes.spines.right": False, "legend.frameon": False})
VT = ["nhit_chan", "sum_energy", "e_over_nhit", "moliere", "zbary", "bar_x", "bar_y"]
VARS = [("nhit_chan", "hits", "core"), ("sum_energy", "energy [MIP]", "core"), ("e_over_nhit", "E/hit", "median"),
        ("moliere", "Moliere radius", "median"), ("zbary", "z barycentre", "median")]
C_VAR = ["#29a3dc", "#c0392b", "#eda100", "#1baf7a", "#4a3aa7"]


def load(paths, acq=False):
    cols = {b: [] for b in VT + (["n_events_acq", "spill", "bcid", "entry"] if acq else [])}
    for p in paths:
        v = uproot.open(f"{p}:ecal").arrays(VT, library="np")
        for b in VT:
            cols[b].append(v[b].astype(float))
        if acq:
            e = uproot.open(f"{p.replace('.valtree.root', '.root')}:ecal").arrays(["n_events_acq", "spill", "bcid"],
                                                                                 library="np")
            assert len(e["spill"]) == len(v["nhit_chan"]), p
            for b in ("n_events_acq", "spill", "bcid"):
                cols[b].append(e[b].astype(float))
            cols["entry"].append(np.arange(len(e["spill"])) / len(e["spill"]))
    c = {b: np.concatenate(x) for b, x in cols.items()}
    k = showerlike(c["nhit_chan"])
    return {b: x[k] for b, x in c.items()}


def stat(x, how):
    x = x[np.isfinite(x)]
    if len(x) < 80:
        return np.nan
    return gauss_core(x)[0] if how == "core" else np.median(x)


def classes(d, s):
    """(label, data mask, sim mask) per class."""
    first = np.zeros(len(d["spill"]), bool)
    order = np.lexsort((d["bcid"], d["spill"]))
    sp = d["spill"][order]
    first[order[np.r_[True, sp[1:] != sp[:-1]]]] = True
    cx, cy = np.median(d["bar_x"]), np.median(d["bar_y"])
    rd = np.hypot(d["bar_x"] - cx, d["bar_y"] - cy)
    rs = np.hypot(s["bar_x"] - np.median(s["bar_x"]), s["bar_y"] - np.median(s["bar_y"]))
    allS = np.ones(len(s["bar_x"]), bool)
    n = d["n_events_acq"]
    out = [("all", np.ones(len(n), bool), allS),
           ("alone in acquisition", n == 1, allS), ("2 in acquisition", n == 2, allS),
           ("3-5 in acquisition", (n >= 3) & (n <= 5), allS), (">5 in acquisition", n > 5, allS),
           ("first of a shared acq.", first & (n > 1), allS), ("later in a shared acq.", ~first & (n > 1), allS)]
    for q in range(4):
        out.append((f"run quarter {q + 1}", (d["entry"] >= q / 4) & (d["entry"] < (q + 1) / 4), allS))
    for lo, hi in ((0, 3), (3, 6), (6, 10), (10, 99)):
        out.append((f"r from centre {lo}-{hi} mm", (rd >= lo) & (rd < hi), (rs >= lo) & (rs < hi)))
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("outdir")
    ap.add_argument("--group", action="append", required=True)
    a = ap.parse_args()
    os.makedirs(a.outdir, exist_ok=True)
    groups = {g["name"]: g for g in load_groups(os.environ["FINAL_GROUPS"])}
    fig, axes = plt.subplots(1, len(a.group), figsize=(8.5 * len(a.group), 8.5), squeeze=False)
    lines = []
    for ax, name in zip(axes[0], a.group):
        g = groups[name]
        d = load([f"{RECO}/{r}/ecal_{r}.valtree.root" for r in g["runs"]], acq=True)
        s = load([sim_valtree(g, "digi")])
        cl = classes(d, s)
        lines.append(f"== {name}: {len(d['nhit_chan'])} data / {len(s['nhit_chan'])} sim shower-like events; "
                     "sim/data per class (Gaussian core for hits and energy, median otherwise)")
        lines.append(f"{'class':<26} {'N data':>7} {'frac':>6} " + " ".join(f"{v[1][:13]:>13}" for v in VARS))
        y = np.arange(len(cl))[::-1]
        for j, (v, lab, how) in enumerate(VARS):
            r = []
            for _, md, ms in (c[:3] for c in cl):
                r.append(stat(s[v][ms], how) / stat(d[v][md], how))
            ax.plot(r, y, "o", color=C_VAR[j], ms=7, label=lab)
            for i, x in enumerate(r):
                cl[i] = cl[i] + (x,)
        for c in cl:
            lines.append(f"{c[0]:<26} {c[1].sum():7d} {c[1].mean():6.3f} " + " ".join(f"{x:13.3f}" for x in c[3:]))
        ax.set_yticks(y); ax.set_yticklabels([f"{c[0]} ({100 * c[1].mean():.0f} %)" for c in cl], fontsize=9.5)
        ax.axvline(1, color=INK, lw=0.8)
        ax.set_xlabel("simulation / data"); ax.set_title(name.replace("_", " "), loc="left")
        ax.legend(fontsize=9, loc="lower right")
    fig.tight_layout(); fig.savefig(f"{a.outdir}/event_classes.png", dpi=120); plt.close(fig)
    txt = "\n".join(lines)
    print(txt)
    open(f"{a.outdir}/event_classes.txt", "w").write(txt + "\n")


if __name__ == "__main__":
    main()

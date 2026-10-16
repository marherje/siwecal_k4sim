#!/usr/bin/env python3
"""Beam points of the uniformity scan U2-U50 (th220, 74 GeV e+, runs 73-141) for analysis/beam_points.txt.

Beam centre = median shower barycentre (bar_x, bar_y) of the shower-like events (final_panels.showerlike) of the
reconstructed runs of each position, as for the other points; sigma = the 74 GeV median of beam_spread.txt (the value
every 74 GeV point of beam_points.txt carries). Runs of one position are pooled when their centres agree within --tol
mm; a run that sits elsewhere becomes its own point named <set>_<label>_e<E>_<x>_<y> (beam_points.txt naming rule).
Runs with fewer than --min-events shower-like events are skipped (run 81: 10 events, run 137: 5).

    python3 analysis/uniformity_points.py            # print the lines
    python3 analysis/uniformity_points.py --write    # ... and append the missing ones to beam_points.txt
"""
import argparse
import os
import sys

import numpy as np
import uproot

_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(_HERE))
from analysis.final_panels import RECO, load_groups, showerlike  # noqa: E402

# Position label -> runs (TB2026-Run_Summary.xlsx)
U = {"U2": [73], "U3": [74], "U4": [75], "U5": [76], "U6": [77, 78], "U7": [79], "U8": [80], "U9": [81, 91],
     "U10": [92, 93], "U11": [94], "U12": [95], "U13": [96], "U14": [97], "U15": [98], "U16": [99], "U17": [100],
     "U18": [101], "U33": [119, 120, 121], "U48": [137], "U49": [141], "U50": [138]}
U.update({f"U{n}": [105 + n - 19] for n in range(19, 33)})
U.update({f"U{n}": [122 + n - 34] for n in range(34, 48)})
SET, E, SIGMA = "th220", "74", ("10.7", "5.9")


def centre(run):
    p = f"{RECO}/{run}/ecal_{run}.valtree.root"
    if not os.path.exists(p):
        return None, 0
    a = uproot.open(p)["ecal"].arrays(["nhit_chan", "bar_x", "bar_y"], library="np")
    k = showerlike(a["nhit_chan"]) if len(a["nhit_chan"]) else np.zeros(0, bool)
    return (np.median(a["bar_x"][k]), np.median(a["bar_y"][k])) if k.sum() else None, int(k.sum())


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--points", default=os.path.join(_HERE, "beam_points.txt"))
    ap.add_argument("--tol", type=float, default=10.0)
    ap.add_argument("--min-events", type=int, default=500)
    ap.add_argument("--write", action="store_true")
    a = ap.parse_args()
    have = {g["name"] for g in load_groups(a.points)}
    lines, missing = [], []
    for label in sorted(U, key=lambda s: int(s[1:])):
        cs = []
        for r in U[label]:
            run = f"TB2026CERN_run_{r:06d}"
            c, n = centre(run)
            if c is None and n == 0 and not os.path.exists(f"{RECO}/{run}/ecal_{run}.valtree.root"):
                missing.append(run)
            print(f"# {label:<4} run {r:3d}: {n:7d} shower-like" + (f"  centre ({c[0]:6.1f}, {c[1]:6.1f})" if c else ""))
            if c is not None and n >= a.min_events:
                cs.append((r, c, n))
        if not cs:
            continue
        # main position = the run with most events; runs within tol pooled with it
        ref = max(cs, key=lambda t: t[2])[1]
        main = [t for t in cs if np.hypot(t[1][0] - ref[0], t[1][1] - ref[1]) <= a.tol]
        rest = [t for t in cs if t not in main]
        w = np.array([t[2] for t in main], float)
        x = np.average([t[1][0] for t in main], weights=w); y = np.average([t[1][1] for t in main], weights=w)
        groups = [(f"{SET}_{label}_e{E}", x, y, main)]
        groups += [(f"{SET}_{label}_e{E}_{t[1][0]:.0f}_{t[1][1]:.0f}", t[1][0], t[1][1], [t]) for t in rest]
        for name, x, y, ts in groups:
            runs = ",".join(str(t[0]) for t in ts)
            lines.append(f"{name:<31} {SET}  {label:<6} {E:<6} {x:4.0f}  {y:4.0f}   {SIGMA[0]:>4}   {SIGMA[1]:>4}  {runs}"
                         + ("" if name not in have else "   # already in beam_points"))
    print("\n".join(lines))
    if missing:
        print(f"# NOT RECONSTRUCTED YET: {' '.join(missing)}")
    if a.write:
        if missing:
            raise SystemExit("some runs are not reconstructed yet; not writing")
        new = [l for l in lines if "already in" not in l]
        with open(a.points, "a") as f:
            f.write("# uniformity scan U2-U50 (analysis/uniformity_points.py)\n" + "\n".join(new) + "\n")
        print(f"# appended {len(new)} points to {a.points}")


if __name__ == "__main__":
    main()

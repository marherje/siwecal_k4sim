#!/usr/bin/env python
"""
Shower barycentre of one group: data against the final_v5 simulation (one beam sigma for every energy) and the
final_v6 one (beam sigma measured in data per energy, beam_spread.py), digitised, shower-like events.

    python analysis/beam_spot_check.py <outdir> --v5-group th210_e20_26_-52 --v6-group th210_S1_e20

bar_r is the barycentre's distance to the detector centre, so it follows the beam spot directly; bar_x and bar_y
show which direction the spot was too narrow or too wide in.
"""
import argparse
import os
import sys

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(_HERE))
from analysis.final_panels import COLOR, RECO, draw_hist, load_groups, showerlike  # noqa: E402
import uproot  # noqa: E402

B = "/eos/experiment/drdcalo/siw-ecal/TB2026-06/Simulation/Processed/adc_vs_tb"
VARS = [("bar_x", "Barycentre x [mm]"), ("bar_y", "Barycentre y [mm]"), ("bar_r", "Barycentre radius [mm]")]


def read(paths):
    cols = {b: [] for b in ["nhit_chan"] + [v for v, _ in VARS]}
    for p in paths:
        for a in uproot.iterate(f"{p}:ecal", list(cols), step_size="400 MB", library="np"):
            for b in cols:
                cols[b].append(a[b])
    c = {b: np.concatenate(v).astype(float) for b, v in cols.items()}
    k = showerlike(c["nhit_chan"])
    return {b: v[k] for b, v in c.items()}


def width(v):
    q25, q75 = np.percentile(v, [25, 75])
    return (q75 - q25) / 1.349


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("outdir")
    ap.add_argument("--v5-group", default="th210_e20_26_-52")
    ap.add_argument("--v6-group", default="th210_S1_e20")
    a = ap.parse_args(argv)
    os.makedirs(a.outdir, exist_ok=True)
    g5 = next(x for x in load_groups(os.path.join(_HERE, "final_v5_groups.txt")) if x["name"] == a.v5_group)
    g6 = next(x for x in load_groups(os.path.join(_HERE, "final_v6_groups.txt")) if x["name"] == a.v6_group)
    s = {"data": read([f"{RECO}/{r}/ecal_{r}.valtree.root" for r in g5["runs"]])}
    for tag, g, d in (("v5", g5, "final_v5"), ("v6", g6, "final_v6")):
        vt = [f for f in os.listdir(f"{B}/{d}/pid/{g['name']}_digi") if f.endswith(".valtree.root")]
        s[tag] = read([f"{B}/{d}/pid/{g['name']}_digi/{vt[0]}"])
    style = {"data": "data", "v5": "nodigi", "v6": "digi"}   # draw_hist styles: points, dotted, hatched
    label = {"data": "Data", "v5": "Sim final_v5 (beam 13.75 x 8.25 mm)", "v6": "Sim final_v6 (beam from data)"}
    fig, axes = plt.subplots(1, 3, figsize=(19, 5.6))
    lines = []
    for ax, (v, xl) in zip(axes, VARS):
        allx = np.concatenate([s[k][v] for k in s])
        lo, hi = np.percentile(allx, [1, 99])
        bins = np.linspace(lo, hi, 41)
        for k in ("v5", "v6", "data"):
            x = s[k][v]
            stat = f"median {np.median(x):.1f}, width {width(x):.1f} mm"
            draw_hist(ax, x, bins, style[k], f"{label[k]}: {stat}")
            lines.append(f"{v} {k} N={len(x)} median={np.median(x):.2f} width={width(x):.2f}")
        ax.set_xlabel(xl, fontsize=12)
        ax.set_ylabel("Events (unit area)", fontsize=11)
        ax.grid(alpha=0.25)
        ax.legend(fontsize=8.5, loc="upper left")
        ax.set_ylim(top=ax.get_ylim()[1] * 1.3)
    runs = ", ".join(("eudaq " if "eudaq" in r else "run ") + str(int(r.split("_")[-1])) for r in g5["runs"])
    fig.suptitle(f"{g5['E']} GeV e$^\\pm$, {g5['set']}, data: {runs}; digitised simulations; shower-like events; "
                 f"N data / v5 / v6 = {len(s['data']['bar_r'])} / {len(s['v5']['bar_r'])} / {len(s['v6']['bar_r'])}",
                 fontsize=12)
    fig.tight_layout(rect=[0, 0, 1, 0.94])
    out = f"{a.outdir}/beam_spot_{a.v6_group}.png"
    fig.savefig(out, dpi=120)
    open(f"{a.outdir}/beam_spot_{a.v6_group}.txt", "w").write("\n".join(lines) + "\n")
    print(out)
    print("\n".join(lines))


if __name__ == "__main__":
    main()

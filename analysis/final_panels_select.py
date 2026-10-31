#!/usr/bin/env python
"""
A reduced 2 x 2 shower panel for one group of the final TB2026 sample: number of hits, calibrated energy, Molière
radius and longitudinal barycentre; data against the simulation without and with digitisation (final_v7).

    python analysis/final_panels_select.py <outdir> [--group th210_e20_26_-52] [--nbins 28]
        [--vars nhit_chan,gev,e_over_nhit,moliere,zbary,mip_likeness]   (any of final_panels.PANELS; 2 x 3 for six)
        [--cal-file <final_panels outdir>/mip_gev_line.txt] [--border-cut MM]

Same inputs, selection and MIP -> GeV line as final_panels.py (valtrees, 0.5 MIP hit cut, shower-like events), with
coarser bins (about half of final_panels.py's 55) so that the 1000-event simulation does not jump bin to bin; for
the integer number of hits the bin width is a whole number of hits, so no bin collects one more integer than its
neighbours. The MIP -> GeV line comes from --cal-file (final_panels.py's mip_gev_line.txt, fitted without the
border cut); without it, the final_v7 line of the report (CAL below). --border-cut MM as in final_panels.py.
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
from analysis.final_panels import (RECO, KINDS, LABEL, COLOR, load_groups, load_valtrees, sim_valtree,  # noqa: E402
                                   gauss_core, core_curve, draw_hist, PANELS as ALL_PANELS, FIT_PANELS, read_cal,
                                   selection_label)

# MIP -> GeV line per sample kind, th210 scan 7.5-74 GeV (final_panels_summary.txt of the report, v3.3; final_v7)
CAL = {"data": (81.92, 42.1), "nodigi": (88.75, 0.0), "digi": (87.15, 2.2)}
PANELS = [("nhit_chan", "Number of hits", True), ("gev", "Calibrated energy [GeV]", True),
          ("moliere", "Molière radius (90 %) [mm]", False), ("zbary", "Longitudinal barycentre [layer]", False)]


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("outdir")
    ap.add_argument("--group", default="th210_e20_26_-52")
    ap.add_argument("--nbins", type=int, default=28)
    ap.add_argument("--vars", default=None, help="comma-separated variables (default: the four of the 2 x 2 panel); "
                    "Gaussian-core fit for the energy-like ones (final_panels.FIT_PANELS), median for the rest")
    ap.add_argument("--cal-file", default=None, help="mip_gev_line.txt of a final_panels.py pass without the cut")
    ap.add_argument("--border-cut", type=float, default=None,
                    help="keep barxy_min_border_distance > MM among the shower-like events")
    a = ap.parse_args(argv)
    if a.border_cut is not None and not a.cal_file:
        ap.error("--border-cut needs --cal-file: the MIP -> GeV line is fitted without the cut")
    os.makedirs(a.outdir, exist_ok=True)
    g = next(x for x in load_groups() if x["name"] == a.group)
    cal = read_cal(a.cal_file) if a.cal_file else CAL
    s = {"data": load_valtrees([f"{RECO}/{r}/ecal_{r}.valtree.root" for r in g["runs"]], a.border_cut)}
    for k in ("nodigi", "digi"):
        s[k] = load_valtrees([sim_valtree(g, k)], a.border_cut)
    for k in KINDS:
        s[k]["gev"] = (s[k]["sum_w"] - cal[k][1]) / cal[k][0]

    panels = PANELS
    if a.vars:
        names = dict(ALL_PANELS)
        panels = [(v, names[v], v in FIT_PANELS) for v in a.vars.split(",")]
    ncol = 2 if len(panels) <= 4 else 3
    nrow = -(-len(panels) // ncol)
    fig, axes = plt.subplots(nrow, ncol, figsize=(6.5 * ncol, 4.8 * nrow))
    lines = []
    for ax, (v, xl, fit) in zip(axes.ravel(), panels):
        xs = {k: s[k][v][np.isfinite(s[k][v])] for k in KINDS}
        # range from each sample's 3-98 % span: the data's tails (low energy, wide Molière: a few % at 74 GeV) must
        # not eat the core's bins; histograms are unit area within this range
        lo = min(np.percentile(x, 3) for x in xs.values())
        hi = max(np.percentile(x, 98) for x in xs.values())
        if v == "nhit_chan":
            w = max(1, int(round((hi - lo) / a.nbins)))
            bins = np.arange(np.floor(lo) - 0.5, np.ceil(hi) + w, w)
        else:
            bins = np.linspace(lo, hi, a.nbins + 1)
        for k in KINDS:
            x = xs[k]
            if fit:
                mu, sg = gauss_core(x)
                fmt = ".1f" if v == "gev" else ".0f"
                draw_hist(ax, x, bins, k, f"{LABEL[k]}: $\\mu$ = {mu:{fmt}}, $\\sigma/\\mu$ = {100 * sg / mu:.1f} %")
                xx, yy = core_curve(x, bins, mu, sg)
                ax.plot(xx, yy, color=COLOR[k], lw=1.3, ls="-" if k != "nodigi" else "--")
                lines.append(f"{v} {k} mu={mu:.4g} res={sg / mu:.4f} N={len(x)}")
            else:
                med = float(np.median(x))
                draw_hist(ax, x, bins, k, f"{LABEL[k]} (median {med:.2f})")
                lines.append(f"{v} {k} median={med:.4g} N={len(x)}")
        if v == "gev":
            ax.axvline(float(g["E"]), color="0.45", lw=1.2, ls="--", label=f"beam energy {g['E']} GeV")
        ax.set_xlabel(xl, fontsize=12)
        ax.set_ylabel("Events (unit area)", fontsize=11)
        ax.grid(alpha=0.25)
        ax.legend(fontsize=8.6, loc="upper left", title="Gaussian core fit" if fit else None, title_fontsize=8.6)
        ax.set_ylim(top=ax.get_ylim()[1] * 1.25)
    nruns = ", ".join(("eudaq " if "eudaq" in r else "run ") + str(int(r.split("_")[-1])) for r in g["runs"])
    fig.suptitle(f"{g['E']} GeV e$^\\pm$, {g['set']}, beam at ({g['pos'].replace('_', ', ')}) mm  —  data: {nruns}; "
                 f"{selection_label(a.border_cut)}, N data / no-digi / digi = {len(xs['data'])} / {len(xs['nodigi'])} / "
                 f"{len(xs['digi'])}", fontsize=12.5)
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    stem = f"panel{len(panels)}_{g['name']}"
    out = f"{a.outdir}/{stem}.png"
    fig.savefig(out, dpi=120)
    open(f"{a.outdir}/{stem}.txt", "w").write("\n".join(lines) + "\n")
    print(out); print("\n".join(lines))


if __name__ == "__main__":
    main()

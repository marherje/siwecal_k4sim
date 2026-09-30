#!/usr/bin/env python
"""
Per-run checks of a reconstructed campaign (default Reconstructed_final) and the beam centroid of each run.

    python analysis/verify_final_campaign.py [--reco DIR] [--out table.txt] [--runs r1,r2,...]

For every run: events, fraction alone in their acquisition (single_run), mean n_events_acq, fraction of hits
with hit_bit == 1, and for shower-like events (nhit > 0.5 x p90, ecal tree) the Gaussian-core mu and sigma/mu of
the energy [MIP] for ALL events and for single_run == 1 only (after the SCA fix the two must agree), plus the
median shower barycentre (bar_x, bar_y) from the PID file -- the beam position of the run's sim sample.
"""
import argparse
import glob
import os
import sys

import numpy as np
import uproot

_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(_HERE))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(_HERE)), "siwecal-tb2026"))
from analysis.compare_event_level import gauss_core  # noqa: E402

B = "/eos/experiment/drdcalo/siw-ecal/TB2026-06"


def run_stats(rdir, run):
    ecal = f"{rdir}/{run}/ecal_{run}.root"
    br = ["nhit_chan", "sum_energy", "single_run", "n_events_acq"]
    cols = {b: [] for b in br}
    nbit = ntot = 0
    for a in uproot.iterate(f"{ecal}:ecal", br + ["hit_bit"], step_size="300 MB", library="np"):
        for b in br:
            cols[b].append(a[b])
        hb = np.concatenate(a["hit_bit"]) if len(a["hit_bit"]) else np.array([])
        nbit += int((hb == 1).sum()); ntot += len(hb)
    c = {b: np.concatenate(v).astype(float) for b, v in cols.items()}
    sh = c["nhit_chan"] > 0.5 * np.percentile(c["nhit_chan"], 90)
    mu_a, sg_a = gauss_core(c["sum_energy"][sh])
    alone = sh & (c["single_run"] == 1)
    mu_s, sg_s = gauss_core(c["sum_energy"][alone]) if alone.sum() >= 50 else (np.nan, np.nan)
    from siwecal_common.edm4hep_pid import PidFileReader
    p = PidFileReader(f"{rdir}/{run}/ecal_{run}.edm4hep.root").scalar_columns()
    nh = np.asarray(p["nhit"], float)
    ps = nh > 0.5 * np.percentile(nh, 90)
    fd = np.asarray(p["fractal_dimension"], float)[ps] if "fractal_dimension" in p else np.array([np.nan])
    return {"N": len(c["nhit_chan"]), "Nsh": int(sh.sum()), "alone": float(np.mean(c["single_run"] == 1)),
            "nacq": float(np.mean(c["n_events_acq"])), "bit": nbit / max(ntot, 1),
            "mu": mu_a, "res": sg_a / mu_a, "mu_s": mu_s, "res_s": sg_s / mu_s if mu_s == mu_s else np.nan,
            "bx": float(np.median(np.asarray(p["bar_x"], float)[ps])),
            "by": float(np.median(np.asarray(p["bar_y"], float)[ps])),
            "nhit": float(np.median(nh[ps])), "fd": float(np.nanmedian(fd))}


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--reco", default=f"{B}/Reconstructed_final")
    ap.add_argument("--runs", default=None)
    ap.add_argument("--out", default=None)
    a = ap.parse_args(argv)
    runs = a.runs.split(",") if a.runs else sorted(os.path.basename(d) for d in glob.glob(f"{a.reco}/TB2026CERN_*"))
    head = (f"{'run':<22}{'events':>9}{'shower':>8}{'alone':>7}{'<nacq>':>7}{'bit=1':>7}{'mu[MIP]':>9}{'s/mu':>7}"
            f"{'mu alone':>9}{'s/mu al':>8}{'bar_x':>7}{'bar_y':>7}{'nhit':>6}{'FD':>6}")
    lines = [f"Campaign {a.reco}; shower-like = nhit > 0.5 x p90; energy = sum_energy [MIP], Gaussian core", head]
    print("\n".join(lines), flush=True)
    for run in runs:
        try:
            s = run_stats(a.reco, run)
        except Exception as e:  # keep going; report the run
            lines.append(f"{run:<22} ERROR {e}"); print(lines[-1], flush=True); continue
        lines.append(f"{run.replace('TB2026CERN_', ''):<22}{s['N']:9d}{s['Nsh']:8d}{s['alone']:7.2f}{s['nacq']:7.2f}"
                     f"{s['bit']:7.3f}{s['mu']:9.0f}{s['res']:7.3f}{s['mu_s']:9.0f}{s['res_s']:8.3f}"
                     f"{s['bx']:7.1f}{s['by']:7.1f}{s['nhit']:6.0f}{s['fd']:6.2f}")
        print(lines[-1], flush=True)
    if a.out:
        open(a.out, "w").write("\n".join(lines) + "\n")


if __name__ == "__main__":
    main()

#!/usr/bin/env python
"""
Does the slow recovery of the slab baseline after earlier events of the acquisition explain the layer term of the sag?
Linear fits of the baseline B of signal-free channels on cm_sag_collect.py records (all runs together):

  base      k Q_chip + b Q_slab                      (same-time charge only)
  +prev     ... + c Q_slab_acq_prev                  (slab charge of the earlier columns of the acquisition)
  +acq      ... + c Q_slab_acq                       (slab charge of the whole acquisition, earlier and later)
  +prev,n   ... + c Q_slab_acq_prev + e n_prev_slab
  each also with the chip x slab interaction d Q_chip Q_slab.

Judged by chi2/ndf of the mean residual in (Q_chip, Q_slab) bins per run, and by the offsets left per (run, slab) when
those are fitted freely on top of the model: if the acquisition terms carry the layer term, the offsets shrink to the
real pedestal drifts.

    python analysis/sag_acq_fit.py <dir with cm_sag_*.npz> <outdir>
"""
import argparse
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from analysis.sag_model_fit import load, judge, RUNINFO  # noqa: E402


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("indir")
    ap.add_argument("outdir")
    a = ap.parse_args()
    os.makedirs(a.outdir, exist_ok=True)
    c = load(a.indir)
    for k in ("q_slab_acq_prev", "q_slab_acq"):
        c[k] = np.maximum(c[k], 0.0)
    B, qc, qs = c["base_far"], c["q_chip"], c["q_slab"]
    qp, qa, npv = c["q_slab_acq_prev"], c["q_slab_acq"], c["n_prev_slab"]
    runs = [r for r in RUNINFO if (c["run"] == r).any()]
    grp = np.array([runs.index(r) for r in c["run"]]) * 15 + c["slab"].astype(int)
    ng = len(runs) * 15
    cnt = np.bincount(grp, minlength=ng)
    models = {"base": [qc, qs], "+prev": [qc, qs, qp], "+acq": [qc, qs, qa], "+prev,n": [qc, qs, qp, npv]}
    for name in list(models):
        models[name + " +QcQs"] = models[name] + [qc * qs * 1e-4]
    lines = [f"{len(B)} records, runs {runs}",
             f"{'model':18s}{'total':>7s}" + "".join(f"{r:>7s}" for r in runs) + "  | free (run,slab) offsets: total chi2,"
             " median |a|, max |a| | parameters [1e-3]"]
    for name, cols in models.items():
        X = np.stack(cols, 1)
        p = np.linalg.lstsq(X, B, rcond=None)[0]
        per, tot = judge(c, X @ p)
        # same model with a free offset per (run, slab): profile the offsets out by demeaning
        dm = lambda x: x - (np.bincount(grp, weights=x, minlength=ng) / np.maximum(cnt, 1))[grp]
        p2 = np.linalg.lstsq(np.stack([dm(x) for x in cols], 1), dm(B), rcond=None)[0]
        off = np.bincount(grp, weights=B - X @ p2, minlength=ng) / np.maximum(cnt, 1)
        _, tot2 = judge(c, off[grp] + X @ p2)
        used = cnt > 1000
        lines.append(f"{name:18s}{tot:7.0f}" + "".join(f"{per[r][0]:7.0f}" for r in runs)
                     + f"  | {tot2:6.0f} {np.median(np.abs(off[used])):5.2f} {np.abs(off[used]).max():5.2f} | "
                     + " ".join(f"{x * 1e3:.4g}" for x in p))
    txt = "\n".join(lines)
    print(txt)
    open(f"{a.outdir}/sag_acq_fit.txt", "w").write(txt + "\n")


if __name__ == "__main__":
    main()

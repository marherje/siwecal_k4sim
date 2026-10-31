#!/usr/bin/env python
"""
Fit the chip baseline sag B (signal-free channels of a chip-SCA) with B = k_chip * Q_chip + f(Q_slab) [+ g(Q_det)],
on the records of cm_sag_collect.py, all runs together, and judge each model by whether it describes every run.

Models (fitted on all records of all runs; k per chip means one k for each (slab, chip id)):
  M0  k Q_chip                          the digitiser now (one k)
  M1  a + k Q_chip                      the pedestal measurement
  M2  k Q_chip + b Q_slab               a linear slab term
  M3  k_chip Q_chip + f(Q_slab)         per-chip k, f = b Q_slab | b ln(1 + Q_slab/q0) | b Q_slab^p
  M4  M3 (best f) + c Q_det             does the rest of the detector add anything
Per-chip k and the slab term are fitted by alternating least squares (k_chip | f, then f | k_chip).
Figure of merit: the record residuals averaged in bins of (Q_chip, Q_slab), per run; chi2/ndf with the statistical
error of each bin mean, per run and in total, and the mean residual per run.

    python analysis/sag_model_fit.py <dir with cm_sag_*.npz> <outdir>

Writes sag_model_fit.txt, sag_model_fit.npz (parameters, k per chip) and sag_model_fit.png.
"""
import argparse
import glob
import os
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from analysis import plot_style as ps  # noqa: E402

RUNINFO = {"e291": "th210 20 GeV", "e287": "th210 52 GeV", "e286": "th210 74 GeV", "20": "th230 20 GeV",
           "13": "th230 52 GeV", "12": "th230 74 GeV", "67": "th220 20 GeV", "72": "th220 74 GeV"}
QC_EDGES = np.r_[0, np.logspace(np.log10(30), np.log10(15000), 13)]
QS_EDGES = np.r_[0, np.logspace(np.log10(300), np.log10(40000), 10)]
QMAX = 15000.0


def load(d):
    recs, runs = [], []
    for p in sorted(glob.glob(f"{d}/cm_sag_*.npz")):
        run = os.path.basename(p)[7:-4]
        z = np.load(p)
        cols = list(z["rec_cols"])
        r = z["rec"]
        keep = np.isfinite(r[:, cols.index("base_far")]) & (r[:, cols.index("q_chip")] < QMAX)
        recs.append(r[keep]); runs += [run] * int(keep.sum())
    r = np.concatenate(recs)
    c = {k: r[:, cols.index(k)].astype(np.float64) for k in cols}
    for k in ("q_chip", "q_slab", "q_det"):    # a fired channel can have a negative low-gain residual
        c[k] = np.maximum(c[k], 0.0)
    c["run"] = np.array(runs)
    c["chipkey"] = (c["slab"].astype(int) * 16 + c["chip"].astype(int))
    return c


MIN_CHIP_RECORDS = 2000   # chips with fewer records (or no charge above 1000 ADC) take the common k


def fit_per_chip_k(c, y):
    """k per chip, no intercept: y ~ k_chip Q_chip; chips with too few records take the common k."""
    key, q = c["chipkey"], c["q_chip"]
    num = np.bincount(key, weights=q * y, minlength=240)
    den = np.bincount(key, weights=q * q, minlength=240)
    n = np.bincount(key, minlength=240)
    nhi = np.bincount(key, weights=(q > 1000).astype(float), minlength=240)
    common = num.sum() / den.sum()
    ok = (n >= MIN_CHIP_RECORDS) & (nhi >= 50)
    return np.where(ok, num / np.maximum(den, 1e-12), common)


def slab_basis(kind, qs, par):
    if kind == "lin":
        return qs
    if kind == "log":
        return np.log1p(qs / par)
    return np.power(np.maximum(qs, 0), par)          # pow


def fit_model(c, name, kind=None, use_det=False):
    """Returns (prediction, parameter dict)."""
    B, qc, qs, qd = c["base_far"], c["q_chip"], c["q_slab"], c["q_det"]
    if name == "M0":
        k = (qc * B).sum() / (qc * qc).sum()
        return k * qc, {"k": k}
    if name == "M1":
        X = np.c_[np.ones_like(qc), qc]
        p = np.linalg.lstsq(X, B, rcond=None)[0]
        return X @ p, {"a": p[0], "k": p[1]}
    if name == "M2":
        X = np.c_[qc, qs]
        p = np.linalg.lstsq(X, B, rcond=None)[0]
        return X @ p, {"k": p[0], "b": p[1]}
    # M3 / M4: alternating least squares, nonlinear parameter of f by a grid
    grid = {"lin": [None], "log": [300, 1000, 2000, 4000, 8000, 16000], "pow": [0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 1.0]}[kind]
    best = None
    for par in grid:
        f = slab_basis(kind, qs, par)
        cols = [f] + ([qd] if use_det else [])
        kchip = np.zeros(240)
        for _ in range(6):
            yk = B - kchip[c["chipkey"]] * qc
            X = np.c_[tuple(cols)]
            p = np.linalg.lstsq(X, yk, rcond=None)[0]
            kchip = fit_per_chip_k(c, B - X @ p)
        pred = kchip[c["chipkey"]] * qc + X @ p
        rss = ((B - pred) ** 2).sum()
        if best is None or rss < best[0]:
            best = (rss, pred, {"f": kind, "par": par, "b": p[0], "c_det": p[1] if use_det else 0.0, "kchip": kchip})
    return best[1], best[2]


def judge(c, pred):
    """chi2/ndf of binned mean residuals per run and overall, and mean residual per run."""
    res = c["base_far"] - pred
    ic = np.clip(np.digitize(c["q_chip"], QC_EDGES) - 1, 0, len(QC_EDGES) - 2)
    js = np.clip(np.digitize(c["q_slab"], QS_EDGES) - 1, 0, len(QS_EDGES) - 2)
    out, chi_all, n_all = {}, 0.0, 0
    for run in RUNINFO:
        m = c["run"] == run
        if not m.any():
            continue
        b = ic[m] * (len(QS_EDGES) - 1) + js[m]
        nb = (len(QC_EDGES) - 1) * (len(QS_EDGES) - 1)
        n = np.bincount(b, minlength=nb); s = np.bincount(b, weights=res[m], minlength=nb)
        s2 = np.bincount(b, weights=res[m] ** 2, minlength=nb)
        ok = n >= 30
        mean = s[ok] / n[ok]; var = s2[ok] / n[ok] - mean ** 2
        chi = (mean ** 2 / (var / n[ok])).sum()
        out[run] = (chi / ok.sum(), res[m].mean())
        chi_all += chi; n_all += ok.sum()
    return out, chi_all / n_all


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("indir")
    ap.add_argument("outdir")
    a = ap.parse_args()
    os.makedirs(a.outdir, exist_ok=True)
    c = load(a.indir)
    runs = [r for r in RUNINFO if (c["run"] == r).any()]
    print(f"{len(c['base_far'])} records, runs {runs}")
    models = [("M0", None, False), ("M1", None, False), ("M2", None, False), ("M3", "lin", False),
              ("M3", "log", False), ("M3", "pow", False)]
    results, lines = {}, []
    for name, kind, det in models:
        pred, par = fit_model(c, name, kind, det)
        per, tot = judge(c, pred)
        results[(name, kind)] = (pred, par, per, tot)
    bestkind = min(("lin", "log", "pow"), key=lambda k: results[("M3", k)][3])
    pred, par = fit_model(c, "M4", bestkind, True)
    results[("M4", bestkind)] = (pred, par, *judge(c, pred))
    head = f"{'model':28s}{'total':>8s}" + "".join(f"{r:>9s}" for r in runs)
    lines += ["chi2/ndf of the mean residual in (Q_chip, Q_slab) bins, per run (statistical errors only)", head]
    for (name, kind), (pred, par, per, tot) in results.items():
        lab = name + (f" f={kind}" if kind else "")
        lines.append(f"{lab:28s}{tot:8.1f}" + "".join(f"{per[r][0]:9.1f}" for r in runs))
    lines += ["", "mean residual per run [ADC] (offset left over)", head.replace("total", "     ")]
    for (name, kind), (pred, par, per, tot) in results.items():
        lab = name + (f" f={kind}" if kind else "")
        lines.append(f"{lab:28s}{'':8s}" + "".join(f"{per[r][1]:9.3f}" for r in runs))
    lines += ["", "parameters:"]
    for (name, kind), (pred, par, per, tot) in results.items():
        shown = {k: (f"median {np.median(v[v != 0]) * 1e3:.3f}e-3, 16-84% {np.percentile(v[v != 0], 16) * 1e3:.3f} .. "
                     f"{np.percentile(v[v != 0], 84) * 1e3:.3f}e-3" if k == "kchip" else f"{v:.4g}" if isinstance(v, float) else v)
                 for k, v in par.items()}
        lines.append(f"  {name} {kind or ''}: {shown}")
    txt = "\n".join(lines)
    print(txt)
    open(f"{a.outdir}/sag_model_fit.txt", "w").write(txt + "\n")
    best = min(results, key=lambda k: results[k][3])
    np.savez(f"{a.outdir}/sag_model_fit.npz", best=str(best), **{f"{n}_{k}__{p}": np.atleast_1d(v)
             for (n, k), (_, par, _, _) in results.items() for p, v in par.items() if not isinstance(v, str)})

    # figure: mean B and model against Q_chip in Q_slab slices, for one run per threshold set, M0 vs the best model
    plt.rcParams.update(ps.RC)
    show = [r for r in ("e286", "12", "72") if r in runs] or runs[:3]
    fig, ax = plt.subplots(1, len(show), figsize=(5.2 * len(show), 4.4), squeeze=False)
    slices = [(0, 2000), (2000, 8000), (8000, 1e9)]
    blues = (ps.DATA, "#1f78a8", "#7cc7ea")
    for j, run in enumerate(show):
        m = c["run"] == run
        for (lo, hi), col, mk in zip(slices, blues, ("o", "s", "^")):
            s = m & (c["q_slab"] >= lo) & (c["q_slab"] < hi)
            ic = np.digitize(c["q_chip"][s], QC_EDGES) - 1
            for lab, pred, sty in (("data", None, None), ("M0 k Q_chip", results[("M0", None)][0], ps.reds(2)[0]),
                                   (f"{best[0]} {best[1]}", results[best][0], ps.reds(2)[1])):
                vals = c["base_far"][s] if pred is None else pred[s]
                n = np.bincount(ic, minlength=len(QC_EDGES)); sm = np.bincount(ic, weights=vals, minlength=len(QC_EDGES))
                qm = np.bincount(ic, weights=c["q_chip"][s], minlength=len(QC_EDGES))
                ok = n >= 30
                x, y = qm[ok] / n[ok], sm[ok] / n[ok]
                tag = f"Q_slab {lo / 1e3:.0f}-{hi / 1e3:.0f}k" if hi < 1e9 else f"Q_slab > {lo / 1e3:.0f}k"
                if pred is None:
                    s2 = np.bincount(ic, weights=vals ** 2, minlength=len(QC_EDGES))
                    err = np.sqrt(np.maximum(s2[ok] / n[ok] - y ** 2, 0) / n[ok])
                    ps.data_points(ax[0, j], x, y, err, label=f"data, {tag}", color=col, fmt=mk)
                else:
                    ax[0, j].plot(x, y, color=sty["color"], ls=sty["ls"], lw=1.6,
                                  label=lab if (lo, hi) == slices[0] else None)
        ax[0, j].set_xscale("log"); ax[0, j].set_title(f"run {run}, {RUNINFO[run]}", fontsize=10)
        ax[0, j].set_xlabel("Q_chip [ADC]")
    ax[0, 0].set_ylabel("<B> of signal-free channels [ADC]")
    ax[0, 0].legend(fontsize=7)
    fig.tight_layout(); fig.savefig(f"{a.outdir}/sag_model_fit.png", dpi=120)


if __name__ == "__main__":
    main()

#!/usr/bin/env python
"""
Measured baseline-sag map B(Q_chip, Q_slab) per threshold set, for AdcDigitizer: instead of a formula, the mean
baseline of signal-free channels in bins of the chip-SCA charge and of the charge of the rest of the slab, from the
cm_sag_collect.py records of the electron runs of that set. The digitiser interpolates bilinearly in
(ln(1 + Q_chip), ln(1 + Q_slab)) between bin centres and clamps at the edges (no extrapolation).

    python analysis/sag_table.py <dir with cm_sag_*.npz> <outdir> [--check]

Writes occupancy_sag_<set>.txt (format below) and, with --check, the per-run residuals and sag_table.png.
File format: comment lines '#'; 'qchip <centres>'; 'qslab <centres>'; then one row of B [ADC] per Q_chip centre.
Bins with fewer than MIN_N records are filled from the nearest filled bin of the same Q_chip row (then column).
"""
import argparse
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from analysis.sag_model_fit import load, judge  # noqa: E402

SETS = {"th210": ["e291", "e287", "e286"], "th220": ["67", "72"], "th230": ["20", "13", "12"]}
QC_EDGES = np.r_[0, np.logspace(np.log10(20), np.log10(15000), 16)]
QS_EDGES = np.r_[0, np.logspace(np.log10(200), np.log10(40000), 12)]
MIN_N = 50


def build(qc, qs, B):
    i = np.clip(np.digitize(qc, QC_EDGES) - 1, 0, len(QC_EDGES) - 2)
    j = np.clip(np.digitize(qs, QS_EDGES) - 1, 0, len(QS_EDGES) - 2)
    shape = (len(QC_EDGES) - 1, len(QS_EDGES) - 1)
    n = np.zeros(shape); s = np.zeros(shape); cx = np.zeros(shape); cy = np.zeros(shape)
    np.add.at(n, (i, j), 1); np.add.at(s, (i, j), B); np.add.at(cx, (i, j), qc); np.add.at(cy, (i, j), qs)
    ok = n >= MIN_N
    tab = np.where(ok, s / np.maximum(n, 1), np.nan)
    # bin centres: mean charge of the records, per axis (over the filled bins)
    xc = np.array([cx[r][ok[r]].sum() / n[r][ok[r]].sum() if ok[r].any() else np.nan for r in range(shape[0])])
    yc = np.array([cy[:, k][ok[:, k]].sum() / n[:, k][ok[:, k]].sum() if ok[:, k].any() else np.nan
                   for k in range(shape[1])])
    keep_r, keep_c = ~np.isnan(xc), ~np.isnan(yc)
    tab, xc, yc, ok = tab[keep_r][:, keep_c], xc[keep_r], yc[keep_c], ok[keep_r][:, keep_c]
    for r in range(tab.shape[0]):                  # fill empty bins along Q_slab from the nearest filled one
        f = np.nonzero(ok[r])[0]
        if len(f):
            for k in np.nonzero(~ok[r])[0]:
                tab[r, k] = tab[r, f[np.argmin(np.abs(f - k))]]
    for k in range(tab.shape[1]):                  # rows still empty: along Q_chip
        f = np.nonzero(~np.isnan(tab[:, k]))[0]
        for r in np.nonzero(np.isnan(tab[:, k]))[0]:
            tab[r, k] = tab[f[np.argmin(np.abs(f - r))], k]
    return xc, yc, tab


def interp(xc, yc, tab, qc, qs):
    """Bilinear in (ln(1+Q_chip), ln(1+Q_slab)), clamped at the edges: what AdcDigitizer does."""
    lx, ly = np.log1p(xc), np.log1p(yc)
    u = np.clip(np.log1p(qc), lx[0], lx[-1]); v = np.clip(np.log1p(qs), ly[0], ly[-1])
    i = np.clip(np.searchsorted(lx, u) - 1, 0, len(lx) - 2); j = np.clip(np.searchsorted(ly, v) - 1, 0, len(ly) - 2)
    tu = (u - lx[i]) / (lx[i + 1] - lx[i]); tv = (v - ly[j]) / (ly[j + 1] - ly[j])
    return ((1 - tu) * (1 - tv) * tab[i, j] + tu * (1 - tv) * tab[i + 1, j]
            + (1 - tu) * tv * tab[i, j + 1] + tu * tv * tab[i + 1, j + 1])


def write(path, th, runs, xc, yc, tab):
    with open(path, "w") as f:
        f.write(f"# Baseline sag B [HG ADC] of signal-free channels, {th}, from electron runs {', '.join(runs)}\n"
                "# (analysis/sag_table.py on cm_sag_collect.py records). Rows: Q_chip centres; columns: Q_slab centres\n"
                "# (charge of the other chips of the slab at the same BCID +-1). Interpolate bilinearly in ln(1+Q),\n"
                "# clamp at the edges.\n")
        f.write("qchip " + " ".join(f"{x:.1f}" for x in xc) + "\n")
        f.write("qslab " + " ".join(f"{y:.1f}" for y in yc) + "\n")
        for row in tab:
            f.write(" ".join(f"{b:.3f}" for b in row) + "\n")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("indir")
    ap.add_argument("outdir")
    ap.add_argument("--check", action="store_true")
    a = ap.parse_args()
    os.makedirs(a.outdir, exist_ok=True)
    c = load(a.indir)
    tables = {}
    for th, runs in SETS.items():
        m = np.isin(c["run"], runs)
        xc, yc, tab = build(c["q_chip"][m], c["q_slab"][m], c["base_far"][m])
        write(f"{a.outdir}/occupancy_sag_{th}.txt", th, runs, xc, yc, tab)
        tables[th] = (xc, yc, tab)
        print(f"{th}: table {tab.shape[0]} x {tab.shape[1]}, B from {tab.max():.2f} to {tab.min():.2f} ADC")
        if a.check:
            cc = {k: v[m] for k, v in c.items()}
            pred = interp(xc, yc, tab, cc["q_chip"], cc["q_slab"])
            per, tot = judge(cc, pred)
            rms = np.sqrt(np.mean((cc["base_far"] - pred) ** 2))
            print(f"   chi2/ndf {tot:.0f} (per run " + ", ".join(f"{r} {per[r][0]:.0f}" for r in runs)
                  + "); mean residual per run " + " ".join(f"{r} {per[r][1]:+.2f}" for r in runs))
    np.save(f"{a.outdir}/sag_tables.npy", tables, allow_pickle=True)


if __name__ == "__main__":
    main()

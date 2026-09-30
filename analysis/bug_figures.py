#!/usr/bin/env python
"""
Figures of the SCA column-order bug for the final artifact (English labels).

    python analysis/bug_figures.py <outdir>

bug_impact.png   per-run energy resolution and energy scale of shower-like events, legacy decoding (all events and
                 single-run events only) vs fixed decoding (all events), from analysis/verify_final_campaign.py
                 tables: final_v5/verify_legacy_hitbit_th210.txt and Reconstructed_final/verify_final.txt.
eudaq_check.png  retrigger signature in the EUDAQ LCIO decoding of run 291, as written by EUDAQ and after the
                 column reversal, against the raw-binary decodings (final_v5/eudaq_check/eudaq_check.txt numbers).
"""
import os
import re
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

B = "/eos/experiment/drdcalo/siw-ecal/TB2026-06"
V5 = f"{B}/Simulation/Processed/adc_vs_tb/final_v5"
C_LEG, C_LEG_ALONE, C_FIX = "#9aa5ad", "#eda100", "#29a3dc"
ENERGY = {"12": 74, "13": 52, "15": 34, "20": 20, "21": 10, "22": 7.5, "41": 7.5, "42": 10, "43": 20, "44": 34,
          "45": 52, "46": 74, "58": 99, "67": 20, "68": 74, "70": 74, "71": 74, "72": 74, "e285": 99, "e286": 74,
          "e287": 52, "e288": 34, "e291": 20, "e292": 10, "e293": 7.5}
SET = {**{k: "th230" for k in ("12", "13", "15", "20", "21", "22", "41", "42", "43", "44", "45", "46", "58")},
       **{k: "th220" for k in ("67", "68", "70", "71", "72")},
       **{k: "th210" for k in ("e285", "e286", "e287", "e288", "e291", "e292", "e293")}}


def read_verify(path):
    rows = {}
    for line in open(path):
        p = line.split()
        if not p or not re.match(r"(eudaq_)?run_\d+", p[0]):
            continue
        num = str(int(p[0].split("_")[-1]))
        key = ("e" if p[0].startswith("eudaq") else "") + num
        f = lambda s: float(s) if s != "nan" else np.nan  # noqa: E731
        rows[key] = {"mu": f(p[6]), "res": f(p[7]), "mu_s": f(p[8]), "res_s": f(p[9])}
    return rows


def bug_impact(out):
    leg = read_verify(f"{V5}/verify_legacy_hitbit_th210.txt")
    fix = read_verify(f"{B}/Reconstructed_final/verify_final.txt")
    keys = [k for k in ENERGY if k in leg and k in fix]
    fig, ax = plt.subplots(1, 2, figsize=(15, 5.6))
    mark = {"th210": "o", "th220": "s", "th230": "^"}
    jit = {"th210": -0.03, "th220": 0.0, "th230": 0.03}
    for k in keys:
        E = ENERGY[k] * (1 + jit[SET[k]]); m = mark[SET[k]]
        ax[0].plot(E, 100 * leg[k]["res"], m, color=C_LEG, ms=8)
        ax[0].plot(E, 100 * leg[k]["res_s"], m, color=C_LEG_ALONE, ms=8, mfc="none", mew=1.8)
        ax[0].plot(E, 100 * fix[k]["res"], m, color=C_FIX, ms=8)
        ref = fix[k]["mu"]
        ax[1].plot(E, leg[k]["mu"] / ref, m, color=C_LEG, ms=8)
        ax[1].plot(E, leg[k]["mu_s"] / ref, m, color=C_LEG_ALONE, ms=8, mfc="none", mew=1.8)
    from matplotlib.lines import Line2D
    h = [Line2D([], [], color=C_LEG, marker="o", ls="none", ms=8, label="legacy decoding, all events"),
         Line2D([], [], color=C_LEG_ALONE, marker="o", ls="none", ms=8, mfc="none", mew=1.8,
                label="legacy decoding, events alone in their acquisition"),
         Line2D([], [], color=C_FIX, marker="o", ls="none", ms=8, label="fixed decoding, all events")]
    h += [Line2D([], [], color="0.35", marker=mark[s], ls="none", ms=8, label=s) for s in ("th210", "th220", "th230")]
    for a in ax:
        a.set_xscale("log"); a.set_xticks([7.5, 10, 20, 34, 52, 74, 99]); a.set_xticklabels(["7.5", "10", "20", "34", "52", "74", "99"])
        a.set_xlabel("Beam energy [GeV]"); a.grid(alpha=0.25, which="both")
    ax[0].set_ylabel("Energy resolution σ/μ [%] (Gaussian core, energy in MIP)")
    ax[1].set_ylabel("Energy scale μ, relative to the fixed decoding")
    ax[1].axhline(1, color="0.5", ls="--", lw=1)
    ax[0].legend(handles=h, fontsize=9, loc="upper right")
    ax[0].set_title("Resolution per run: the width came from the pairing, not the detector", fontsize=11)
    ax[1].set_title("Energy scale per run: legacy all-event mean is 5-20 % low", fontsize=11)
    fig.suptitle("Same runs, same event building and calibration family, shower-like e± events", fontsize=13)
    fig.tight_layout(rect=[0, 0, 1, 0.94])
    fig.savefig(f"{out}/bug_impact.png", dpi=110); plt.close(fig)


def eudaq_check(out):
    # chunk_copy_check.py on 3000 acquisitions of run 291 (EUDAQ LCIO) and chunk 0 of the raw decodings
    cats = ["n=2\ngap 1-3", "n=2\nfar", "n=3\ngap 1-3", "n=3\nfar", "n=4\ngap 1-3", "n=4\nfar"]
    as_written = [0.189, 0.000, 0.005, 0.086, 0.100, 0.092]
    corrected = [0.189, 0.000, 0.238, 0.001, 0.252, 0.002]
    raw_legacy = [0.155, 0.000, 0.004, 0.106, 0.063, 0.087]
    raw_fixed = [0.155, 0.000, 0.245, 0.002, 0.145, 0.002]
    x = np.arange(len(cats)); w = 0.2
    fig, ax = plt.subplots(figsize=(12, 5.2))
    for i, (v, c, lab, hatch) in enumerate(((raw_legacy, C_LEG, "raw binary, legacy decoder", None),
                                            (as_written, "#c0392b", "EUDAQ LCIO as written by EUDAQ", None),
                                            (raw_fixed, C_FIX, "raw binary, fixed decoder", None),
                                            (corrected, "#1f5f86", "EUDAQ LCIO after column reversal", "//"))):
        ax.bar(x + (i - 1.5) * w, v, w * 0.92, color=c, label=lab, hatch=hatch, edgecolor="white" if hatch else None)
    ax.set_xticks(x); ax.set_xticklabels(cats)
    ax.set_ylabel("Fraction of adjacent column pairs\nwith the retrigger signature")
    ax.set_xlabel("Columns in the chip (n) and BCID distance of the pair (gap 1-3 = consecutive trigger + retrigger)")
    ax.legend(fontsize=9); ax.grid(axis="y", alpha=0.25)
    ax.set_title("Run 291: the retrigger copy must sit in BCID-consecutive pairs. EUDAQ's LCIO puts it in far pairs,\n"
                 "exactly like the legacy raw decoder (800/800 acquisitions bit-identical); reversed, it matches the fix",
                 fontsize=10.5)
    fig.tight_layout(); fig.savefig(f"{out}/eudaq_check.png", dpi=110); plt.close(fig)


if __name__ == "__main__":
    out = sys.argv[1]
    os.makedirs(out, exist_ok=True)
    bug_impact(out)
    eudaq_check(out)

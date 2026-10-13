#!/usr/bin/env python
"""
The chip baseline drop chip by chip, and whether it is large enough to explain the energy-per-hit loss in data.

Input: the per chip-SCA records of analysis/cm_trigger_test.py (cm_trigger_<run>.npz, SCA columns that follow another
one of the same chip within 1-3 BCIDs vetoed: there the channels fired in the earlier column still carry their signal
without a hit bit and pull the baseline up): hit bits, charge Q of the fired channels (low gain / k, HG ADC) and
baseline B (median ADC - pedestal of the non-fired channels >= 2 pads from any fired one, HG ADC).

  cm_chip_<slab>_<chip>.png  B against Q and against the hit bits for one chip, 20 and 74 GeV, with the all-chip median
  cm_slope_map.png           slope dB/dQ [ADC per 1000 ADC] of every chip (pooled e291, e288, e287, e286, e166), map
                             per slab, and its spread
  cm_magnitude.png           for listed chips: the fractional energy-per-hit loss the baseline drop implies,
                             (B - B_quiet) / <Q per fired channel>, against hit bits, next to the measured <E/hit>
                             data/sim - 1 (analysis/chip_occupancy.py output)

    python analysis/cm_per_chip.py <cm_trigger npz dir> <outdir> [--occ <chip_occupancy.txt>]
"""
import argparse
import os
import re

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

RUNS = {"e291": "20 GeV", "e288": "34 GeV", "e287": "52 GeV", "e286": "74 GeV", "e166": "99 GeV"}
NB = np.array([1, 3, 6, 10, 15, 20, 25, 30, 40, 65])
Q_EDGES = np.array([0, 100, 200, 350, 500, 750, 1000, 1500, 2000, 3000, 4000, 6000, 8000, 12000, 20000])
B_QUIET = -0.6        # muon runs, 0 hit bits: what the pedestal tables already contain
CHIPS = [(6, 11), (6, 10), (5, 11), (7, 11), (4, 11)]
C_CHIP = ["#c0392b", "#eda100", "#29a3dc", "#1baf7a", "#4a3aa7"]
INK, MUTED, GRID = "#1b2227", "#56626b", "#d7dde2"
plt.rcParams.update({"font.size": 11, "axes.edgecolor": MUTED, "axes.labelcolor": INK, "xtick.color": MUTED,
                     "ytick.color": MUTED, "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.7,
                     "axes.spines.top": False, "axes.spines.right": False, "legend.frameon": False})


def load(d, r):
    """Per chip-SCA records of analysis/cm_trigger_test.py (retrigger columns vetoed); slab 6 chip 11 from its own,
    larger sample when there is one."""
    z = np.load(f"{d}/cm_trigger_{r}.npz")
    R = {c: z["rec"][:, i] for i, c in enumerate(z["rec_cols"])}
    if os.path.exists(f"{d}/cm_trigger_{r}_6_11.npz"):
        z = np.load(f"{d}/cm_trigger_{r}_6_11.npz")
        S = {c: z["rec"][:, i] for i, c in enumerate(z["rec_cols"])}
        keep = ~((R["slab"] == 6) & (R["chip"] == 11))
        R = {k: np.concatenate([R[k][keep], S[k]]) for k in R}
    ok = np.isfinite(R["base_far"]) & (R["n_far"] >= 3)
    return {k: v[ok] for k, v in R.items()}


def binned(x, y, edges, nmin=15):
    k = np.searchsorted(edges, x, side="right") - 1
    med, ctr, err = [], [], []
    for i in range(len(edges) - 1):
        s = k == i
        if s.sum() >= nmin:
            med.append(np.median(y[s])); ctr.append(np.median(x[s]))
            err.append(1.2533 * np.std(y[s]) / np.sqrt(s.sum()))
    return np.array(ctr), np.array(med), np.array(err)


def slope(q, b, nmin=40):
    """dB/dQ from a straight line through the binned medians (Q up to 8000 ADC), ADC per 1000 ADC."""
    x, y, e = binned(q, b, Q_EDGES[Q_EDGES <= 8000], nmin=nmin)
    if len(x) < 4 or x.max() < 1500:
        return np.nan, np.nan
    w = 1 / np.maximum(e, 0.05) ** 2
    A = np.vstack([np.ones_like(x), x / 1000]).T
    cov = np.linalg.inv(A.T @ (A * w[:, None]))
    p = cov @ (A.T @ (w * y))
    return p[1], np.sqrt(cov[1, 1])


def read_occ(path):
    """{(label, slab, chip): E/hit data/sim per NB bin} from chip_occupancy.txt."""
    out, lab = {}, None
    for line in open(path):
        m = re.match(r"## th210 S1 (\d+) GeV", line)
        if m:
            lab = f"{m.group(1)} GeV"
            continue
        m = re.match(r"slab (\d+) chip (\d+): (.*?)\s*\|", line)
        if m and lab:
            out[(lab, int(m.group(1)), int(m.group(2)))] = np.array([float(v) for v in m.group(3).split()])
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("indir")
    ap.add_argument("outdir")
    ap.add_argument("--occ", default=None)
    a = ap.parse_args()
    D = {r: load(a.indir, r) for r in RUNS}
    lines = []

    # 1. one chip against the all-chip median, Q and hit bits
    fig, ax = plt.subplots(1, 2, figsize=(12, 4.6))
    for r, ls in (("e291", "-"), ("e286", "--")):
        R = D[r]
        x, y, _ = binned(R["q"], R["base_far"], Q_EDGES, nmin=50)
        ax[0].plot(x, y, ls, color=MUTED, lw=1.6, label=f"all chips, {RUNS[r]}")
        x, y, _ = binned(R["nbits"], R["base_far"], NB, nmin=50)
        ax[1].plot(x, y, ls, color=MUTED, lw=1.6, label=f"all chips, {RUNS[r]}")
        for (s, c), col in zip(CHIPS[:2], C_CHIP):
            m = (R["slab"] == s) & (R["chip"] == c)
            x, y, e = binned(R["q"][m], R["base_far"][m], Q_EDGES, nmin=12)
            ax[0].errorbar(x, y, e, fmt=ls + "o", ms=4, color=col, label=f"slab {s} chip {c}, {RUNS[r]}")
            x, y, e = binned(R["nbits"][m], R["base_far"][m], NB, nmin=12)
            ax[1].errorbar(x, y, e, fmt=ls + "o", ms=4, color=col)
    ax[0].set_xscale("log"); ax[0].set_xlabel("charge of the fired channels in the chip-SCA [HG ADC]")
    ax[1].set_xlabel("channels with a hit bit in the chip-SCA")
    for x_ in ax:
        x_.set_ylabel("baseline of non-fired channels [ADC]")
        x_.axhline(B_QUIET, color=INK, lw=0.8, ls=":")
    ax[0].legend(fontsize=8.5, ncol=2)
    fig.suptitle("Baseline (ADC - pedestal of channels without a hit bit, >= 2 pads from any hit), median per bin",
                 fontsize=11)
    fig.tight_layout(); fig.savefig(f"{a.outdir}/cm_chip_6_11.png", dpi=130); plt.close(fig)

    # 2. slope map, pooled energies
    P = {k: np.concatenate([D[r][k] for r in ("e291", "e288", "e287", "e286", "e166")]) for k in ("slab", "chip", "q", "base_far")}
    S, E = np.full((15, 16), np.nan), np.full((15, 16), np.nan)
    for s in range(15):
        for c in range(16):
            m = (P["slab"] == s) & (P["chip"] == c)
            if m.sum() > 200:
                S[s, c], E[s, c] = slope(P["q"][m], P["base_far"][m])
    g, ge = slope(P["q"], P["base_far"], nmin=200)
    ok = np.isfinite(S) & (E < 0.3)
    lines.append(f"Slope dB/dQ [ADC per 1000 ADC], pooled 20-99 GeV: all chips {g:+.3f} +- {ge:.3f}; "
                 f"{ok.sum()} chips measured (error < 0.3): median {np.median(S[ok]):+.3f}, "
                 f"16-84 % {np.percentile(S[ok], 16):+.3f} .. {np.percentile(S[ok], 84):+.3f}, "
                 f"median error {np.median(E[ok]):.3f}")
    lines.append(f"  slab 6 chip 11: {S[6, 11]:+.3f} +- {E[6, 11]:.3f}; slab 6 chip 10: {S[6, 10]:+.3f} +- {E[6, 10]:.3f}")
    pull = (S[ok] - np.median(S[ok])) / E[ok]
    lines.append(f"  chip-to-chip spread beyond errors: RMS of pulls {np.sqrt(np.mean(pull ** 2)):.2f} "
                 f"(1 = consistent with one common slope)")
    lines.append("  per slab median: " + " ".join(f"{np.nanmedian(np.where(ok[s], S[s], np.nan)):+.2f}"
                                                  if ok[s].any() else "  .  " for s in range(15)))
    fig, ax = plt.subplots(1, 2, figsize=(12, 4.6), gridspec_kw={"width_ratios": [1.5, 1]})
    Sm = np.where(ok, S, np.nan)
    im = ax[0].imshow(Sm.T, origin="lower", cmap="viridis_r", vmin=-2.0, vmax=0, aspect="auto")
    ax[0].set_xlabel("slab"); ax[0].set_ylabel("chip"); ax[0].set_xticks(range(15)); ax[0].set_yticks(range(0, 16, 3))
    ax[0].grid(False)
    ax[0].add_patch(plt.Rectangle((5.5, 10.5), 1, 1, fill=False, ec="#c0392b", lw=2))
    fig.colorbar(im, ax=ax[0], label="dB/dQ [ADC per 1000 ADC]")
    ax[1].hist(S[ok], bins=np.arange(-2.5, 0.75, 0.1), color="#2a78d6", alpha=0.85)
    ax[1].axvline(g, color=INK, ls="--", lw=1, label=f"all chips {g:+.2f}")
    ax[1].axvline(S[6, 11], color="#c0392b", lw=2, label=f"slab 6 chip 11 {S[6, 11]:+.2f}")
    ax[1].set_xlabel("dB/dQ [ADC per 1000 ADC]"); ax[1].set_ylabel("chips"); ax[1].legend(fontsize=9)
    fig.suptitle("Baseline slope per chip against the charge in the chip-SCA, 20-99 GeV pooled (grey: < 200 chip-SCAs "
                 "or error > 0.3)", fontsize=11)
    fig.tight_layout(); fig.savefig(f"{a.outdir}/cm_slope_map.png", dpi=130); plt.close(fig)

    # 3. magnitude: implied E/hit loss against measured
    occ = read_occ(a.occ) if a.occ else {}
    fig, ax = plt.subplots(1, 2, figsize=(12, 5.4), sharey=True)
    nbc = 0.5 * (NB[1:] + NB[:-1])
    for i, (r, lab) in enumerate((("e291", "20 GeV"), ("e286", "74 GeV"))):
        R = D[r]
        lines.append(f"== {lab}: implied E/hit change (B - B_quiet) / <Q per fired channel> | measured E/hit data/sim - 1")
        for (s, c), col in zip(CHIPS, C_CHIP):
            m = (R["slab"] == s) & (R["chip"] == c)
            kb = np.searchsorted(NB, R["nbits"][m], side="right") - 1
            imp = np.full(len(nbc), np.nan)
            for j in range(len(nbc)):
                t = kb == j
                if t.sum() >= 10:
                    imp[j] = (np.median(R["base_far"][m][t]) - B_QUIET) / np.median(R["q"][m][t] / R["nbits"][m][t])
            ax[i].plot(nbc, 100 * imp, "o-", color=col, label=f"slab {s} chip {c}: baseline drop")
            meas = occ.get((lab, s, c))
            if meas is not None:
                ax[i].plot(nbc, 100 * (meas - 1), "s:", color=col, mfc="none", label="measured E/hit data/sim - 1")
            lines.append(f"  slab {s} chip {c:2d}: " + " ".join(f"{100 * v:+5.1f}" if np.isfinite(v) else "   . "
                                                               for v in imp)
                         + (" | " + " ".join(f"{100 * (v - 1):+5.1f}" if np.isfinite(v) else "   . " for v in meas)
                            if meas is not None else ""))
        ax[i].axhline(0, color=INK, lw=0.8)
        ax[i].set_title(lab); ax[i].set_xlabel("hits (hit bits) in the chip")
    ax[0].set_ylabel("energy per hit change [%]")
    fig.legend(*ax[0].get_legend_handles_labels(), loc="lower center", ncol=5, fontsize=8.5)
    fig.suptitle("Filled: what the measured baseline drop does to the energy per hit. Open: the loss data/sim shows.",
                 fontsize=11)
    fig.tight_layout(rect=(0, 0.1, 1, 1)); fig.savefig(f"{a.outdir}/cm_magnitude.png", dpi=130); plt.close(fig)
    lines.append(f"  bins (hit bits): {' '.join(f'{lo}-{hi - 1}' for lo, hi in zip(NB[:-1], NB[1:]))}")
    txt = "\n".join(lines)
    print(txt)
    open(f"{a.outdir}/cm_per_chip.txt", "w").write(txt + "\n")


if __name__ == "__main__":
    main()

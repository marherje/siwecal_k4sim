#!/usr/bin/env python
"""
Figures of analysis/cm_mip_peak.py against the digitised simulation: position of the 1-MIP peak and low/high gain ratio
against the charge Q of the other fired channels of the chip-SCA.

    python analysis/cm_mip_peak_plot.py <dir> "<label>=<data run>|<sim ecal tree>" [...]

<data run> is the suffix of <dir>/cm_mip_peak_<run>.npz; <sim ecal tree> an ecal tree (hit_slab/chip/energy/hg/lg) or
'-' for none. In the simulation a chip has one trigger per event, so Q is the energy of the other hits of the same
slab and chip. Writes <dir>/cm_mip_peak.png and cm_mip_peak.txt.
"""
import sys

import numpy as np
import uproot
import awkward as ak
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from scipy.optimize import curve_fit  # noqa: E402

INK, MUTED, GRID = "#1b2227", "#56626b", "#d7dde2"
RAMP = ["#c6dbef", "#9ecae1", "#6baed6", "#4292c6", "#2171b5", "#08519c", "#08306b", "#041c40"]
C_RUN = ["#29a3dc", "#eda100", "#c0392b"]
plt.rcParams.update({"font.size": 11, "axes.edgecolor": MUTED, "axes.labelcolor": INK, "xtick.color": MUTED,
                     "ytick.color": MUTED, "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.7,
                     "axes.spines.top": False, "axes.spines.right": False, "legend.frameon": False})


def peak(h, bins, lo=0.6, hi=2.2):
    """Mode of the MIP peak: Gaussian fitted around the maximum in [lo, hi] (window -0.2, +0.3 MIP)."""
    ctr = 0.5 * (bins[1:] + bins[:-1])
    w = (ctr > lo) & (ctr < hi)
    if h[w].sum() < 400:
        return np.nan, np.nan
    hs = np.convolve(h, np.ones(5) / 5, mode="same")
    m0 = ctr[w][np.argmax(hs[w])]
    f = (ctr > m0 - 0.2) & (ctr < m0 + 0.3)
    try:
        p, c = curve_fit(lambda x, a, m, s: a * np.exp(-0.5 * ((x - m) / s) ** 2), ctr[f], h[f],
                         p0=[h[f].max(), m0, 0.3], sigma=np.sqrt(np.maximum(h[f], 1)))
        return (p[1], np.sqrt(c[1, 1])) if lo < p[1] < hi else (np.nan, np.nan)
    except RuntimeError:
        return np.nan, np.nan


def ratio_median(h, bins):
    if h.sum() < 200:
        return np.nan
    c = np.cumsum(h) / h.sum()
    return np.interp(0.5, c, 0.5 * (bins[1:] + bins[:-1]))


def sim_hist(path, q_edges, s_bins, nmax=30000):
    a = uproot.open(path)["ecal"].arrays(["hit_slab", "hit_chip", "hit_energy", "hit_hg", "hit_lg"], entry_stop=nmax,
                                         library="ak")
    key = a["hit_slab"] * 16 + a["hit_chip"]
    e = a["hit_energy"]
    # energy of the other hits of the same slab and chip
    flat_k, flat_e = ak.flatten(key), ak.flatten(e)
    ev = np.repeat(np.arange(len(a)), ak.to_numpy(ak.num(e)))
    g = ev * 240 + ak.to_numpy(flat_k)
    tot = np.bincount(g, weights=ak.to_numpy(flat_e))[g]
    q = tot - ak.to_numpy(flat_e)
    s = ak.to_numpy(flat_e)
    nq = len(q_edges) - 1
    kq = np.clip(np.searchsorted(q_edges, q, side="right") - 1, 0, nq - 1)
    return np.stack([np.histogram(s[kq == k], bins=s_bins)[0] for k in range(nq)])


def main():
    d, items = sys.argv[1], [x.split("=", 1) for x in sys.argv[2:]]
    fig, ax = plt.subplots(2, 2, figsize=(14, 10.5))
    lines = []
    for i, (lab, spec) in enumerate(items):
        run, sim = spec.split("|")
        z = np.load(f"{d}/cm_mip_peak_{run}.npz")
        qe, sb, rb = z["q_edges"], z["s_bins"], z["r_bins"]
        ctr = 0.5 * (sb[1:] + sb[:-1])
        qc = np.sqrt(np.maximum(qe[:-1], 1) * np.minimum(qe[1:], 2 * np.maximum(qe[:-1], 1)))
        hs = sim_hist(sim, qe, sb) if sim != "-" else None
        res = {k: np.array([peak(z[k][j], sb) for j in range(len(qe) - 1)]) for k in ("raw_all", "cor_all")}
        res["sim"] = (np.array([peak(hs[j], sb) for j in range(len(qe) - 1)]) if hs is not None
                      else np.full((len(qe) - 1, 2), np.nan))
        rat = np.array([ratio_median(z["ratio"][j], rb) for j in range(len(qe) - 1)])
        rat_raw = np.array([ratio_median(z["ratio_raw"][j], rb) for j in range(len(qe) - 1)])
        lines.append(f"== {lab} (data {run}, {int(z['n_acq'])} acquisitions, baseline measured in "
                     f"{float(z['b_measured_frac']):.2f} of the chip-SCAs; sim {sim})")
        lines.append(f"{'Q [MIP]':>11} {'N data':>10} {'peak raw':>9} {'peak cor':>9} {'peak sim':>9} "
                     f"{'cor/sim':>8} {'LG/HG':>7} {'raw':>7} {'B_LG':>6}")
        for j in range(len(qe) - 1):
            ql = f"{qe[j]:.0f}-{qe[j + 1]:.0f}" if qe[j + 1] < 1e8 else f">{qe[j]:.0f}"
            lines.append(f"{ql:>11} {int(z['cor_all'][j][ctr > 0.6].sum()):10d} {res['raw_all'][j, 0]:9.3f} "
                         f"{res['cor_all'][j, 0]:9.3f} {res['sim'][j, 0]:9.3f} "
                         f"{res['cor_all'][j, 0] / res['sim'][j, 0]:8.3f} {rat[j]:7.4f} {rat_raw[j]:7.4f} {z['b_lg'][j]:+6.2f}")
        col = C_RUN[i % len(C_RUN)]
        ax[0, 1].errorbar(qc, res["cor_all"][:, 0], yerr=res["cor_all"][:, 1], fmt="o-", color=col, ms=6,
                          label=f"data {lab}, baseline-corrected")
        ax[0, 1].plot(qc, res["raw_all"][:, 0], "o:", color=col, mfc="none", ms=6, label=f"data {lab}, raw")
        if hs is not None:
            ax[0, 1].errorbar(qc, res["sim"][:, 0], yerr=res["sim"][:, 1], fmt="s--", color=col, mfc="white", ms=6,
                              label=f"sim {lab}")
        ax[1, 1].plot(qc, rat, "o-", color=col, ms=6, label=f"data {lab}, LG baseline-corrected")
        ax[1, 1].plot(qc, rat_raw, "o:", color=col, mfc="none", ms=6, label=f"data {lab}, LG raw")
        if i == len(items) - 1:   # spectra of the highest-energy sample
            for j in range(len(qe) - 1):
                w = (ctr > 0.6) & (ctr < 3.5)
                hd = z["cor_all"][j]
                if hd[w].sum() < 400:
                    continue
                ql = f"{qe[j]:.0f}-{qe[j + 1]:.0f}" if qe[j + 1] < 1e8 else f">{qe[j]:.0f}"
                ax[0, 0].plot(ctr[w], hd[w] / hd[w].sum() / np.diff(sb)[0], color=RAMP[j], lw=1.6, label=f"Q {ql}")
                if hs is not None and hs[j][w].sum() > 400:
                    ax[1, 0].plot(ctr[w], hs[j][w] / hs[j][w].sum() / np.diff(sb)[0], color=RAMP[j], lw=1.6,
                                  label=f"Q {ql}")
            ax[0, 0].set_title(f"Data {lab}: all channels, baseline-corrected", loc="left", fontsize=11)
            ax[1, 0].set_title(f"Simulation {lab}: hits", loc="left", fontsize=11)
    for x in (ax[0, 0], ax[1, 0]):
        x.set_xlabel("signal / MIP of the channel (high gain)"); x.set_ylabel("density (0.6-3.5 MIP)")
        x.legend(fontsize=8.5)
    ax[0, 1].set_xscale("log"); ax[1, 1].set_xscale("log")
    ax[0, 1].set_ylabel("1-MIP peak position [MIP]"); ax[0, 1].legend(fontsize=8)
    ax[1, 1].set_ylabel("median (LG - ped - B$_{LG}$) / (HG - ped - B), HG 150-800 ADC"); ax[1, 1].legend(fontsize=8.5)
    for x in (ax[0, 1], ax[1, 1]):
        x.set_xlabel("charge of the other fired channels in the chip-SCA [MIP]")
    fig.tight_layout(); fig.savefig(f"{d}/cm_mip_peak.png", dpi=120); plt.close(fig)
    txt = "\n".join(lines)
    print(txt)
    open(f"{d}/cm_mip_peak.txt", "w").write(txt + "\n")


if __name__ == "__main__":
    main()

#!/usr/bin/env python
"""
The muon gain scan re-read, per threshold set: the Pearson chi2 of muon_gain_scan.py, its expected Poisson noise floor
(n_bins / N_sim_hits), the chi2 with both variances, and an independent observable, the tag-and-probe layer
efficiency of muon tracks (a probe layer is counted when >= 10 other layers fired), data band against the
simulation at each gain. Port of final_v3/plots/gain_scan_floor.py with the inputs as arguments.

    python analysis/gain_scan_floor.py --out gain_scan_floor.png \
        --data th210=<ecal.root> --data th220=... --data th230=... \
        --sims th210=<dir with th210_mu_g*.root> --sims th220=... --sims th230=... \
        [--threshold th210=16.11 ...] [--txt out.txt]
"""
import argparse
import glob
import os
import sys

import numpy as np
import uproot
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from analysis.fit_adc_scale import layers_hit  # noqa: E402

COL = {"th210": "#29a3dc", "th220": "#eda100", "th230": "#c0392b"}
THR = {"th210": 16.11, "th220": 20.23, "th230": 26.07}


def load(path, max_events=None, probe=(4, 5, 6, 7, 8, 9, 10)):
    a = uproot.open(path)["ecal"].arrays(["nhit_chan", "hit_slab", "hit_hg", "hit_ismasked"], entry_stop=max_events,
                                         library="np")
    keep = (a["nhit_chan"] <= 30) & (layers_hit(a["hit_slab"], a["hit_ismasked"]) >= 10)
    hg = np.concatenate(a["hit_hg"][keep]).astype(float)
    m = np.concatenate(a["hit_ismasked"][keep]).astype(int)
    counts, _ = np.histogram(hg[m == 0], bins=np.arange(0, 151, 1))
    num = den = 0
    for n, sl, mm in zip(a["nhit_chan"], a["hit_slab"], a["hit_ismasked"]):
        if n > 30:
            continue
        layers = set(sl[mm == 0].astype(int).tolist())
        for p in probe:
            if len(layers - {p}) >= 10:
                den += 1
                num += (p in layers)
    return counts, num / max(den, 1), np.sqrt(max(num, 1)) / max(den, 1)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", required=True)
    ap.add_argument("--data", action="append", required=True, metavar="TH=ECAL")
    ap.add_argument("--sims", action="append", required=True, metavar="TH=DIR")
    ap.add_argument("--threshold", action="append", default=[], metavar="TH=ADC")
    ap.add_argument("--txt", default=None)
    a = ap.parse_args(argv)
    data = dict(x.split("=", 1) for x in a.data)
    sims = dict(x.split("=", 1) for x in a.sims)
    thr = dict(THR, **{k: float(v) for k, v in (x.split("=", 1) for x in a.threshold)})
    ths = [t for t in ("th210", "th220", "th230") if t in data and t in sims]
    fig, axes = plt.subplots(2, len(ths), figsize=(5.2 * len(ths), 8.4), squeeze=False)
    lines = []
    for col, th in enumerate(ths):
        lo, hi = int(round(thr[th])) + 1, 80
        w = slice(lo, hi)
        cd, ed, eed = load(data[th], 40000)
        Nd = cd[w].sum(); d = cd / Nd
        G, P, F, X, E, EE = [], [], [], [], [], []
        files = sorted(glob.glob(f"{sims[th]}/{th}_mu_g*.root"), key=lambda f: float(f.split("_g")[-1][:-5]))
        for f in files:
            g = float(f.split("_g")[-1][:-5])
            cs, es, ees = load(f)
            Ns = cs[w].sum(); s = cs / Ns
            P.append(np.sum((s[w] - d[w]) ** 2 / np.maximum(d[w], 1e-4)))
            F.append(np.sum(s[w] / Ns / np.maximum(d[w], 1e-4)) + (hi - lo) / Nd)
            X.append(np.sum((s[w] - d[w]) ** 2 / (s[w] / Ns + d[w] / Nd)))
            G.append(g); E.append(es); EE.append(ees)
            lines.append(f"{th} gain {g:5.2f}: chi2(2var) {X[-1]:7.0f}/{hi - lo} bins  probe eff sim {es:.4f} "
                         f"+- {ees:.4f}  data {ed:.4f} +- {eed:.4f}")
        G = np.array(G); E = np.array(E)
        # gain at which the simulated probe efficiency crosses the data (linear interpolation)
        cross = np.nan
        for i in range(len(G) - 1):
            if (E[i] - ed) * (E[i + 1] - ed) <= 0 and E[i + 1] != E[i]:
                cross = G[i] + (ed - E[i]) * (G[i + 1] - G[i]) / (E[i + 1] - E[i]); break
        imin = int(np.argmin(X))
        lines.append(f"{th}: chi2 minimum on the grid at {G[imin]:.2f} ADC/MIP; probe efficiency crosses the data at "
                     f"{cross:.2f} ADC/MIP")
        ax = axes[0, col]
        ax.plot(G, np.array(P) * 1e3, "o-", color=COL[th], lw=1.8, label="Pearson $\\chi^2$ (data as expectation)")
        ax.plot(G, np.array(F) * 1e3, "s--", color="#78828f", lw=1.6, label="its Poisson noise floor")
        ax.set_title(f"{th}: discriminator {thr[th]:.1f} ADC, window [{lo}, {hi}) ADC", fontsize=10.5)
        ax.set_xlabel("adc_per_mip [ADC/MIP]"); ax.set_ylabel("$\\chi^2$ ×10³"); ax.grid(alpha=0.3); ax.legend(fontsize=8)
        ax2 = axes[1, col]
        ax2.plot(G, X, "o-", color=COL[th], lw=1.8, label="$\\chi^2$ with both Poisson variances")
        ax2.axhline(hi - lo, color="grey", ls=":", lw=1.0, label=f"n bins = {hi - lo}")
        ax2.set_ylim(0, max(4 * (hi - lo), min(max(X), 12 * (hi - lo))))
        ax2.set_xlabel("adc_per_mip [ADC/MIP]"); ax2.set_ylabel("$\\chi^2$"); ax2.grid(alpha=0.3)
        ax3 = ax2.twinx()
        ax3.errorbar(G, E, yerr=EE, fmt="^", color="#7a5a00", ms=5, label="sim: probe-layer efficiency")
        ax3.axhspan(ed - eed, ed + eed, color="#7a5a00", alpha=0.18, label=f"data: {ed:.4f}")
        ax3.set_ylabel("tag-and-probe layer efficiency", color="#7a5a00"); ax3.tick_params(axis="y", colors="#7a5a00")
        h1, l1 = ax2.get_legend_handles_labels(); h2, l2 = ax3.get_legend_handles_labels()
        ax2.legend(h1 + h2, l1 + l2, fontsize=8, loc="upper center")
    fig.suptitle("Muon gain scan: shape $\\chi^2$ and tag-and-probe layer efficiency against the simulated gain",
                 fontsize=12)
    fig.tight_layout(rect=[0, 0, 1, 0.96]); fig.savefig(a.out, dpi=110)
    print("\n".join(lines))
    if a.txt:
        open(a.txt, "w").write("\n".join(lines) + "\n")


if __name__ == "__main__":
    main()

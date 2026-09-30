#!/usr/bin/env python
"""
Per-slab trigger model: layer profiles (hits and energy per layer, shower-like events) of the data against the
digitised simulation with the per-slab discriminator overrides (the production final_v5 samples) and without
them (one threshold for every slab, REAL_THRESHOLD_MIP), as ratios sim / data per layer.

    python analysis/perslab_effect.py --out perslab_effect.png \
        --case "th210 20 GeV=<data valtree(s, comma)>|<sim valtree with overrides>|<sim valtree without>" [--case ...]
"""
import argparse

import numpy as np
import uproot
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

C_DATA, C_WITH, C_WITHOUT = "#29a3dc", "#c0392b", "#8c8c8c"


def profiles(paths):
    hits, en = [], []
    for p in paths:
        a = uproot.open(p)["ecal"].arrays(["nhit_chan", "hits_per_layer", "energy_per_layer"], library="np")
        k = a["nhit_chan"] > 0.5 * np.percentile(a["nhit_chan"], 90)
        hits.append(np.stack(a["hits_per_layer"][k])); en.append(np.stack(a["energy_per_layer"][k]))
    return np.concatenate(hits).mean(axis=0), np.concatenate(en).mean(axis=0)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", required=True)
    ap.add_argument("--case", action="append", required=True)
    a = ap.parse_args(argv)
    fig, axes = plt.subplots(2, len(a.case), figsize=(7.2 * len(a.case), 8.4), squeeze=False)
    for col, case in enumerate(a.case):
        label, files = case.split("=", 1)
        data, with_, without = files.split("|")
        hd, ed = profiles(data.split(","))
        hw, ew = profiles([with_])
        ho, eo = profiles([without])
        L = np.arange(len(hd))
        for row, (d, w, o, yl) in enumerate(((hd, hw, ho, "hits per layer, sim / data"),
                                             (ed, ew, eo, "energy per layer [MIP], sim / data"))):
            ax = axes[row, col]
            ax.axhline(1, color="0.5", ls="--", lw=1)
            ax.plot(L, o / d, "s--", color=C_WITHOUT, lw=1.6, ms=6, mfc="white", label="sim, one threshold for all slabs")
            ax.plot(L, w / d, "o-", color=C_WITH, lw=2, ms=7, label="sim, per-slab discriminator (final_v5)")
            ax.axvspan(11.5, 12.5, color=C_DATA, alpha=0.12, lw=0)
            ax.text(12, ax.get_ylim()[1], "slab 12\n(COB, DAC 243)", ha="center", va="top", fontsize=8.5)
            ax.set_xlabel("Layer"); ax.set_ylabel(yl); ax.grid(alpha=0.25)
            if row == 0:
                ax.set_title(label, fontsize=11.5); ax.legend(fontsize=8.5, loc="lower left")
    fig.suptitle("Per-slab trigger model: layer profiles of shower-like events, simulation over data", fontsize=12.5)
    fig.tight_layout(rect=[0, 0, 1, 0.95]); fig.savefig(a.out, dpi=110)


if __name__ == "__main__":
    main()

#!/usr/bin/env python
"""
Does a chip lose energy when many of its channels fire? Mean energy per hit of one chip in an event against the number
of hits that chip has in the event, data over digitised simulation, for a list of (slab, chip) and one or more beam
groups. Shower-like events from the valtrees (0.5 MIP hit cut in both).

    python analysis/chip_occupancy.py --out <dir> --chips 6:11,6:10,5:11,7:11,4:11 \
        --group "20 GeV=<data valtree>|<sim valtree>" --group "74 GeV=<data valtree>|<sim valtree>"

Writes chip_occupancy.png (one panel per group: <E/hit> data/sim against chip hits, one line per chip, plus the
fraction of events in each bin) and chip_occupancy.txt.
"""
import argparse

import numpy as np
import uproot
import awkward as ak
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

from trigger_maps import showerlike, C_GRID

NBINS = np.array([1, 3, 6, 10, 15, 20, 25, 30, 40, 65])
C_CHIP = ["#c0392b", "#eda100", "#29a3dc", "#1baf7a", "#4a3aa7", "#8c8c8c"]


def per_chip(path, chips, nmax=60000):
    a = uproot.open(path)["ecal"].arrays(["nhit_chan", "hit_slab", "hit_chip", "hit_energy"], entry_stop=nmax,
                                         library="ak")
    a = a[showerlike(ak.to_numpy(a["nhit_chan"]).astype(float))]
    out = {}
    for L, c in chips:
        m = (a["hit_slab"] == L) & (a["hit_chip"] == c)
        out[(L, c)] = (ak.to_numpy(ak.sum(m, axis=1)), ak.to_numpy(ak.sum(a["hit_energy"] * m, axis=1)))
    return out


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", required=True)
    ap.add_argument("--chips", required=True)
    ap.add_argument("--group", action="append", required=True)
    a = ap.parse_args(argv)
    chips = [tuple(int(x) for x in s.split(":")) for s in a.chips.split(",")]
    nc = 0.5 * (NBINS[1:] + NBINS[:-1])
    fig, axes = plt.subplots(2, len(a.group), figsize=(8.4 * len(a.group), 9.6), squeeze=False,
                             gridspec_kw={"height_ratios": [3, 1.3]}, constrained_layout=True, sharex=True)
    lines = []
    for col, g in enumerate(a.group):
        label, files = g.split("=", 1)
        dp, sp = files.split("|")
        D, S = per_chip(dp, chips), per_chip(sp, chips)
        lines.append(f"## {label}: <E/hit> data/sim per chip-hit bin {list(NBINS)}")
        for i, ch in enumerate(chips):
            (nd, ed), (ns, es) = D[ch], S[ch]
            r, e, fd, fs = [], [], [], []
            for lo, hi in zip(NBINS[:-1], NBINS[1:]):
                kd = (nd >= lo) & (nd < hi); ks = (ns >= lo) & (ns < hi)
                fd.append(kd.mean()); fs.append(ks.mean())
                if kd.sum() > 30 and ks.sum() > 30:
                    md = ed[kd].sum() / nd[kd].sum(); ms = es[ks].sum() / ns[ks].sum()
                    # error of the mean energy per hit from the event-to-event spread
                    sd = np.std(ed[kd] / nd[kd]) / np.sqrt(kd.sum()); ss = np.std(es[ks] / ns[ks]) / np.sqrt(ks.sum())
                    r.append(md / ms); e.append(md / ms * np.hypot(sd / md, ss / ms))
                else:
                    r.append(np.nan); e.append(np.nan)
            ax = axes[0, col]
            ax.errorbar(nc + 0.25 * (i - len(chips) / 2), r, yerr=e, fmt="o-", color=C_CHIP[i % len(C_CHIP)], lw=2,
                        ms=7, capsize=2, label=f"slab {ch[0]} chip {ch[1]} (data {nd.mean():.1f}, sim {ns.mean():.1f} hits)")
            fd_, fs_ = np.array(fd), np.array(fs)
            ok = (fd_ > 0.002) & (fs_ > 0.002)
            axes[1, col].plot(nc[ok], fd_[ok] / fs_[ok], "o-", color=C_CHIP[i % len(C_CHIP)], lw=1.6, ms=6)
            lines.append(f"slab {ch[0]} chip {ch[1]}: " + " ".join(f"{x:.3f}" for x in r)
                         + "   | event fraction data/sim " + " ".join(f"{x:.2f}" for x in np.array(fd) / np.maximum(np.array(fs), 1e-9)))
        ax = axes[0, col]
        ax.axhline(1, color="0.5", lw=1, ls="--")
        ax.set_ylabel("<energy per hit> of the chip, data / sim", fontsize=13)
        ax.set_title(label, fontsize=14); ax.grid(color=C_GRID, lw=0.6); ax.legend(fontsize=10, frameon=False)
        axes[1, col].axhline(1, color="0.5", lw=1, ls="--")
        axes[1, col].set_yscale("log")
        axes[1, col].set_ylabel("events in bin, data / sim\n(bins with > 0.2 % of events)", fontsize=12)
        axes[1, col].set_xlabel("hits in the chip in the event", fontsize=13); axes[1, col].grid(color=C_GRID, lw=0.6)
    fig.savefig(f"{a.out}/chip_occupancy.png", dpi=110); plt.close(fig)
    with open(f"{a.out}/chip_occupancy.txt", "w") as f:
        f.write("\n".join(lines) + "\n")
    print("\n".join(lines))


if __name__ == "__main__":
    main()

#!/usr/bin/env python
"""
Where does the simulation have more hits than the data? Hits (and energy) per shower-like event in bins of layer and
of the distance r to the event's barycentre axis, digitised and undigitised simulation over data, for a few groups.

    python analysis/hit_deficit.py <outdir> <label>=<data valtree>|<sim digi valtree>|<sim nodigi valtree> [...]

Writes hit_deficit.png (r profile of sim/data per group, and the layer x r map of the digitised ratio for each) and
hit_deficit.txt.
"""
import sys

import numpy as np
import uproot
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

R_EDGES = np.array([0, 5, 10, 15, 20, 30, 40, 60, 90, 200.0])
N_LAYERS = 15


def profile(path, nmax=40000):
    a = uproot.open(path)["ecal"].arrays(["nhit_chan", "bar_x", "bar_y", "hit_x", "hit_y", "hit_slab", "hit_energy"],
                                         entry_stop=nmax, library="np")
    k = a["nhit_chan"] > 0.5 * np.percentile(a["nhit_chan"], 90)
    H = np.zeros((N_LAYERS, len(R_EDGES) - 1)); En = np.zeros_like(H); n = 0
    for bx, by, x, y, s, e in zip(a["bar_x"][k], a["bar_y"][k], a["hit_x"][k], a["hit_y"][k], a["hit_slab"][k],
                                  a["hit_energy"][k]):
        r = np.hypot(x - bx, y - by)
        H += np.histogram2d(s, r, bins=[np.arange(N_LAYERS + 1) - 0.5, R_EDGES])[0]
        En += np.histogram2d(s, r, bins=[np.arange(N_LAYERS + 1) - 0.5, R_EDGES], weights=e)[0]
        n += 1
    return H / n, En / n, n


def main():
    out = sys.argv[1]
    cases = [c.split("=", 1) for c in sys.argv[2:]]
    fig, axes = plt.subplots(2, len(cases), figsize=(6.2 * len(cases), 9.5), squeeze=False)
    rc = 0.5 * (R_EDGES[1:] + R_EDGES[:-1]); rc[-1] = 120
    lines = []
    for col, (label, files) in enumerate(cases):
        dp, sp, np_ = files.split("|")
        Hd, Ed, nd = profile(dp); Hs, Es, ns = profile(sp); Hn, En_, nn = profile(np_)
        ax = axes[0, col]
        for H, Hn_, c, ls, lab in ((Hs, None, "#c0392b", "-", "sim digitised"), (Hn, None, "#eda100", ":", "sim undigitised")):
            ax.plot(rc, H.sum(0) / Hd.sum(0), "o" + ls, color=c, lw=2, ms=6, label=f"{lab}: hits")
        ax.plot(rc, Es.sum(0) / Ed.sum(0), "s--", color="#c0392b", lw=1.4, ms=5, mfc="white", label="sim digitised: energy")
        ax.axhline(1, color="0.5", ls="--", lw=1); ax.set_xscale("log"); ax.set_xticks([2.5, 7.5, 12.5, 17.5, 25, 35, 50, 75, 120])
        ax.set_xticklabels(["2.5", "7.5", "12.5", "17.5", "25", "35", "50", "75", ">90"])
        ax.set_xlabel("distance to the shower axis r [mm]"); ax.set_ylabel("sim / data, per event"); ax.grid(alpha=0.25)
        ax.set_title(label, fontsize=11); ax.legend(fontsize=8.5)
        ax2 = axes[1, col]
        ratio = np.where(Hd > 0.02, Hs / np.maximum(Hd, 1e-9), np.nan)
        im = ax2.imshow(ratio, origin="lower", aspect="auto", cmap="RdBu_r", vmin=0.7, vmax=1.3,
                        extent=[-0.5, len(rc) - 0.5, -0.5, N_LAYERS - 0.5])
        ax2.set_xticks(range(len(rc))); ax2.set_xticklabels([f"{int(a)}-{int(b)}" if b < 200 else f">{int(a)}"
                                                             for a, b in zip(R_EDGES[:-1], R_EDGES[1:])], rotation=45, fontsize=8)
        ax2.set_xlabel("r [mm]"); ax2.set_ylabel("layer"); plt.colorbar(im, ax=ax2, label="hits, sim digitised / data")
        lines.append(f"== {label}: events data {nd}, sim {ns}, nodigi {nn}")
        lines.append("   r bin [mm]      " + "".join(f"{int(a):>7}" for a in R_EDGES[:-1]))
        lines.append("   hits data/ev    " + "".join(f"{v:7.2f}" for v in Hd.sum(0)))
        lines.append("   hits digi/data  " + "".join(f"{v:7.3f}" for v in Hs.sum(0) / Hd.sum(0)))
        lines.append("   hits nodg/data  " + "".join(f"{v:7.3f}" for v in Hn.sum(0) / Hd.sum(0)))
        lines.append("   E digi/data     " + "".join(f"{v:7.3f}" for v in Es.sum(0) / Ed.sum(0)))
        lines.append("   layer hits digi/data " + " ".join(f"{v:.2f}" for v in Hs.sum(1) / Hd.sum(1)))
    fig.suptitle("Hits per shower-like event, simulation over data, against the distance to the shower axis", fontsize=12.5)
    fig.tight_layout(rect=[0, 0, 1, 0.96]); fig.savefig(f"{out}/hit_deficit.png", dpi=105)
    open(f"{out}/hit_deficit.txt", "w").write("\n".join(lines) + "\n"); print("\n".join(lines))


if __name__ == "__main__":
    main()

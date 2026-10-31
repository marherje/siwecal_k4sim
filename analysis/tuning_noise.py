#!/usr/bin/env python
"""
sim_tuning 2.2: noisier / quieter channels and chips, and whole-chip dropouts, data against digitised simulation.

  isolated hits   hits with no other hit within one pad (8 neighbours) in the same layer and > 50 mm from the event
                  barycentre: noise or stray signal, not shower core. Rate per event per cell and per chip, and
                  their energy spectrum. Hot cells: data rate > 5 x simulation and > 1e-3 per event.
  chip dropouts   per chip, the fraction of events in which the chip has no hit, data and simulation, for chips the
                  simulation fills (mean > 3 hits/event); excess dropout probability p = (P_d - P_s) / (1 - P_s).

Writes <out>/noise_<group>.png (maps), dropout_<group>.png, noise.txt and noise_<group>.npz (p_drop[slab, chip],
isolated rates, hot cell list) for the strategies.

    python analysis/tuning_noise.py <outdir> [groups]
"""
import os
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from analysis.sim_tuning import (load, load_groups, data_paths, sim_paths, chip_key, OUT, NL, NG, NCHIP)  # noqa

R_ISO = 50.0


def isolated(s):
    base = (s["ev"].astype(np.int64) * NL + s["slab"]) * (NG + 2)
    key = (base + s["iy"] + 1) * (NG + 2) + s["ix"] + 1
    keys = np.unique(key)
    has_nb = np.zeros(len(key), bool)
    for dy in (-1, 0, 1):
        for dx in (-1, 0, 1):
            if dx or dy:
                has_nb |= np.isin(key + dy * (NG + 2) + dx, keys)
    r = np.hypot(s["x"] - s["bar_x"][s["ev"]], s["y"] - s["bar_y"][s["ev"]])
    return (~has_nb) & (r > R_ISO)


def maps(s):
    iso = isolated(s)
    nev = int(s["nev"])
    cc = (s["slab"].astype(int) * NG + s["iy"]) * NG + s["ix"]
    cell = np.bincount(cc[iso], minlength=NL * NG * NG).reshape(NL, NG, NG) / nev
    ci = s["slab"].astype(int) * NCHIP + s["chip"]
    chip = np.bincount(ci[iso], minlength=NL * NCHIP).reshape(NL, NCHIP) / nev
    # dropouts
    u = np.unique(chip_key(s))
    filled = np.zeros((NL, NCHIP))
    np.add.at(filled, ((u // NCHIP) % NL, u % NCHIP), 1)
    p_empty = 1 - filled / nev
    mean_hits = np.bincount(ci, minlength=NL * NCHIP).reshape(NL, NCHIP) / nev
    return {"iso": iso, "cell": cell, "chip": chip, "p_empty": p_empty, "mean_hits": mean_hits,
            "iso_e": s["e"][iso], "iso_frac": iso.mean(), "iso_per_ev": iso.sum() / nev}


def draw_slabs(fig, gs_row, m, title, vmin, vmax, cmap, nrow=1):
    axs = []
    for l in range(NL):
        a = fig.add_subplot(gs_row[l])
        im = a.imshow(m[l], origin="lower", vmin=vmin, vmax=vmax, cmap=cmap)
        a.set_xticks([]); a.set_yticks([]); a.set_title(f"L{l}", fontsize=7)
        axs.append(a)
    axs[0].set_ylabel(title, fontsize=8)
    return im


def main():
    out = sys.argv[1]
    want = sys.argv[2:] or ["th210_e20_26_-52", "th210_e74_22_-49", "th230_E4_e20", "th230_E4_e74"]
    os.makedirs(out, exist_ok=True)
    groups = {g["name"]: g for g in load_groups()}
    lines = ["group  iso/ev data sim  | iso E median data sim | hot cells | chips with p_drop>2% (slab:chip p)"]
    for n in want:
        g = groups[n]
        d = load(data_paths(g), cache=f"{OUT}/cache/data_{n}.npz")
        s = load(sim_paths(n), cache=f"{OUT}/cache/v7_{n}.npz")
        md, ms = maps(d), maps(s)
        nd, ns = int(d["nev"]), int(s["nev"])
        hot = np.argwhere((md["cell"] > 5 * np.maximum(ms["cell"], 1.0 / ns)) & (md["cell"] > 1e-3))
        pd_, ps = md["p_empty"], ms["p_empty"]
        p_drop = np.where(ms["mean_hits"] > 3, np.clip((pd_ - ps) / np.maximum(1 - ps, 1e-9), 0, 1), np.nan)
        big = np.argwhere(np.nan_to_num(p_drop) > 0.02)
        lines.append(f"{n}  {md['iso_per_ev']:.3f} {ms['iso_per_ev']:.3f} | {np.median(md['iso_e']):.2f} "
                     f"{np.median(ms['iso_e']):.2f} | {len(hot)} | " +
                     " ".join(f"{l}:{c} {p_drop[l, c]:.3f}" for l, c in big))
        print(lines[-1], flush=True)
        np.savez(f"{out}/noise_{n}.npz", p_drop=p_drop, iso_cell_data=md["cell"], iso_cell_sim=ms["cell"],
                 iso_chip_data=md["chip"], iso_chip_sim=ms["chip"], hot=hot, p_empty_data=pd_, p_empty_sim=ps)

        fig = plt.figure(figsize=(18, 6.4))
        gs = fig.add_gridspec(3, NL + 1, width_ratios=[1] * NL + [0.08], hspace=0.35, wspace=0.08)
        lr = np.log10(np.maximum(md["cell"], 1e-5)); ls_ = np.log10(np.maximum(ms["cell"], 1e-5))
        im = draw_slabs(fig, [gs[0, l] for l in range(NL)], lr, "data: log10 iso/ev", -5, -1, "magma")
        fig.colorbar(im, cax=fig.add_subplot(gs[0, NL]))
        im = draw_slabs(fig, [gs[1, l] for l in range(NL)], ls_, "sim: log10 iso/ev", -5, -1, "magma")
        fig.colorbar(im, cax=fig.add_subplot(gs[1, NL]))
        with np.errstate(divide="ignore", invalid="ignore"):
            rr = np.where((md["cell"] * nd >= 5) | (ms["cell"] * ns >= 5),
                          np.log2(np.maximum(md["cell"], 0.5 / nd) / np.maximum(ms["cell"], 0.5 / ns)), np.nan)
        im = draw_slabs(fig, [gs[2, l] for l in range(NL)], rr, "log2 data/sim", -3, 3, "RdBu_r")
        fig.colorbar(im, cax=fig.add_subplot(gs[2, NL]))
        fig.suptitle(f"{n}: isolated hits (no neighbour, > {R_ISO:.0f} mm from the barycentre) per event and cell")
        fig.savefig(f"{out}/noise_{n}.png", dpi=100, bbox_inches="tight"); plt.close(fig)

        fig, ax = plt.subplots(1, 3, figsize=(15, 4.2))
        a = ax[0].imshow(np.where(ms["mean_hits"] > 3, pd_ - ps, np.nan).T, origin="lower", cmap="RdBu_r",
                         vmin=-0.05, vmax=0.05, aspect="auto")
        fig.colorbar(a, ax=ax[0]); ax[0].set_xlabel("layer"); ax[0].set_ylabel("chip")
        ax[0].set_title("P(chip empty): data - sim (chips sim > 3 hits/ev)")
        bins = np.linspace(0.5, 6, 45)
        ax[1].hist(md["iso_e"], bins, density=True, histtype="step", label="data", color="#29a3dc")
        ax[1].hist(ms["iso_e"], bins, density=True, histtype="step", label="sim", color="#c0392b")
        ax[1].set_xlabel("isolated hit energy [MIP]"); ax[1].legend()
        lay = np.arange(NL)
        ax[2].plot(lay, md["chip"].sum(1), "-o", color="#29a3dc", label="data")
        ax[2].plot(lay, ms["chip"].sum(1), "-s", color="#c0392b", label="sim")
        ax[2].set_xlabel("layer"); ax[2].set_ylabel("isolated hits / event"); ax[2].legend()
        fig.suptitle(n); fig.tight_layout(); fig.savefig(f"{out}/dropout_{n}.png", dpi=100); plt.close(fig)
    open(f"{out}/noise.txt", "w").write("\n".join(lines) + "\n")


if __name__ == "__main__":
    main()

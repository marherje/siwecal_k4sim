#!/usr/bin/env python
"""
Where, cell by cell and chip by chip, does the digitised simulation have more hits than the data? Shower-like events
of one beam group, hits per event in every (slab, chip, channel) from the valtrees (the PID stage's 0.5 MIP cut, the
same in data and simulation), for one or more digitisations of the same simulated sample (e.g. the per-slab trigger
model against the per-channel muon trigger tables).

    python analysis/trigger_maps.py --out <dir> --label "th210 S1 20 GeV" \
        --data <data valtree>[,<more>] --sim "per-slab=<valtree>" --sim "muon tables=<valtree>" \
        [--trigger-table <Trigger_channel_highgain.txt>] [--max-events 40000]

Writes, per simulation, into <dir>:
  cells_<tag>_<sim>.png   15 slabs, hits/event sim over data per cell (log2 colour, grey = too few hits)
  chips_<tag>_<sim>.png   15 slabs, the same ratio summed per chip, painted over the chip's cells
and once per group:
  coherence_<tag>.png     per slab: how much of the cell-to-cell scatter of log(sim/data) the chip explains, against
                          the same 8x8 blocks shifted by 4 pads (a smooth, chip-blind pattern explains as much as the
                          chips do) and the statistical floor; log(sim/data) per cell against the cell's muon trigger
                          threshold minus its chip's (does the data lose hits where the discriminator sits high?)
  chipocc_<tag>.png       hits in a chip against the energy the chip collects in the event (data and sims), and
                          the spectrum of the hit energies near the threshold
  trigger_maps_<tag>.txt  the numbers behind them
"""
import argparse
import os

import numpy as np
import uproot
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.colors import TwoSlopeNorm  # noqa: E402

N_SLAB, N_CHIP, N_CHAN = 15, 16, 64
NCELL = N_SLAB * N_CHIP * N_CHAN
C_DATA = "#29a3dc"
C_SIM = ["#c0392b", "#1baf7a", "#eda100", "#4a3aa7"]
C_GRID = "#d8d8d4"
MIN_HITS = 25          # cells (chips) with fewer hits in data or simulation are not drawn
BR = ["nhit_chan", "hit_slab", "hit_chip", "hit_chan", "hit_x", "hit_y", "hit_energy"]


def showerlike(n):
    return n > 0.5 * np.percentile(n, 90)


def load(paths, nmax):
    """Shower-like events: per-cell hit counts, cell positions, and the per-(event, chip) hit count and energy."""
    counts = np.zeros(NCELL)
    xsum, ysum, nxy = np.zeros(NCELL), np.zeros(NCELL), np.zeros(NCELL)
    chip_n, chip_e, hit_e = [], [], []
    nev = 0
    for p in paths:
        for a in uproot.iterate(f"{p}:ecal", BR, step_size="300 MB", library="ak"):
            import awkward as ak
            n = ak.to_numpy(a["nhit_chan"]).astype(float)
            a = a[showerlike(n)]
            if nmax and nev + len(a) > nmax:
                a = a[: nmax - nev]
            m = ak.num(a["hit_slab"])
            cell = (ak.to_numpy(ak.flatten(a["hit_slab"])).astype(int) * N_CHIP
                    + ak.to_numpy(ak.flatten(a["hit_chip"])).astype(int)) * N_CHAN \
                + ak.to_numpy(ak.flatten(a["hit_chan"])).astype(int)
            e = ak.to_numpy(ak.flatten(a["hit_energy"])).astype(float)
            counts += np.bincount(cell, minlength=NCELL)
            xsum += np.bincount(cell, ak.to_numpy(ak.flatten(a["hit_x"])), minlength=NCELL)
            ysum += np.bincount(cell, ak.to_numpy(ak.flatten(a["hit_y"])), minlength=NCELL)
            nxy += np.bincount(cell, minlength=NCELL)
            # one entry per (event, slab, chip) with at least one hit
            ev = np.repeat(np.arange(len(a)), ak.to_numpy(m))
            key = ev * (N_SLAB * N_CHIP) + cell // N_CHAN
            u, inv = np.unique(key, return_inverse=True)
            chip_n.append(np.bincount(inv).astype(float))
            chip_e.append(np.bincount(inv, e))
            hit_e.append(e[e < 4.0])
            nev += len(a)
            if nmax and nev >= nmax:
                break
        if nmax and nev >= nmax:
            break
    with np.errstate(invalid="ignore", divide="ignore"):
        pos = np.stack([xsum / nxy, ysum / nxy], axis=1)
    return {"counts": counts, "nev": nev, "pos": pos, "chip_n": np.concatenate(chip_n),
            "chip_e": np.concatenate(chip_e), "hit_e": np.concatenate(hit_e)}


def read_trigger(path):
    """Trigger_channel_highgain.txt -> mu_adc per cell (NaN where masked) and the fit status."""
    mu = np.full(NCELL, np.nan); st = np.full(NCELL, -1)
    for line in open(path):
        if line.startswith("#"):
            continue
        f = line.split()
        if len(f) < 9:
            continue
        c = (int(f[0]) * N_CHIP + int(f[1])) * N_CHAN + int(f[2])
        st[c] = int(f[8])
        if st[c] != 3:
            mu[c] = float(f[3])
    return mu, st


def pad_index(pos):
    """Pad column/row of every cell from its mean hit position (32 distinct x and y values)."""
    ok = np.isfinite(pos[:, 0])
    ix = np.full(NCELL, -1); iy = np.full(NCELL, -1)
    for k, col in ((0, ix), (1, iy)):
        vals = np.unique(np.round(pos[ok, k], 1))
        col[ok] = np.searchsorted(vals, np.round(pos[ok, k], 1))
    return ix, iy


def ratio_cells(d, s):
    """log2(sim/data) of hits per event per cell, its statistical error, and the cells with enough hits."""
    ok = (d["counts"] >= MIN_HITS) & (s["counts"] >= MIN_HITS)
    with np.errstate(divide="ignore", invalid="ignore"):
        r = (s["counts"] / s["nev"]) / (d["counts"] / d["nev"])
        lr = np.log2(r)
        err = np.sqrt(1 / d["counts"] + 1 / s["counts"]) / np.log(2)
    lr[~ok] = np.nan
    return lr, err, ok


def ratio_chips(d, s):
    dc = d["counts"].reshape(-1, N_CHAN).sum(1); sc = s["counts"].reshape(-1, N_CHAN).sum(1)
    ok = (dc >= MIN_HITS) & (sc >= MIN_HITS)
    with np.errstate(divide="ignore", invalid="ignore"):
        lr = np.log2((sc / s["nev"]) / (dc / d["nev"]))
    lr[~ok] = np.nan
    return lr


def slab_grid(values, ix, iy, slab):
    g = np.full((32, 32), np.nan)
    sl = slice(slab * N_CHIP * N_CHAN, (slab + 1) * N_CHIP * N_CHAN)
    v, x, y = values[sl], ix[sl], iy[sl]
    k = (x >= 0) & np.isfinite(v)
    g[y[k], x[k]] = v[k]
    return g


def chip_outlines(ax, ix, iy, slab):
    """Thin lines between pads of different chips."""
    sl = slice(slab * N_CHIP * N_CHAN, (slab + 1) * N_CHIP * N_CHAN)
    chip = np.full((32, 32), -1)
    x, y = ix[sl], iy[sl]; k = x >= 0
    chip[y[k], x[k]] = (np.arange(N_CHIP * N_CHAN) // N_CHAN)[k]
    for j in range(32):
        for i in range(32):
            if chip[j, i] < 0:
                continue
            if i < 31 and chip[j, i + 1] >= 0 and chip[j, i + 1] != chip[j, i]:
                ax.plot([i + 0.5, i + 0.5], [j - 0.5, j + 0.5], color="0.15", lw=0.6)
            if j < 31 and chip[j + 1, i] >= 0 and chip[j + 1, i] != chip[j, i]:
                ax.plot([i - 0.5, i + 0.5], [j + 0.5, j + 0.5], color="0.15", lw=0.6)


def map_figure(values, ix, iy, title, path, lim=1.0, per_slab_text=None):
    fig, axes = plt.subplots(3, 5, figsize=(22, 13.6), constrained_layout=True)
    norm = TwoSlopeNorm(vcenter=0.0, vmin=-lim, vmax=lim)
    cmap = plt.get_cmap("RdBu_r").copy(); cmap.set_bad("#bdbdbd")
    for slab, ax in enumerate(axes.ravel()):
        g = slab_grid(values, ix, iy, slab)
        im = ax.imshow(g, origin="lower", cmap=cmap, norm=norm, interpolation="nearest")
        chip_outlines(ax, ix, iy, slab)
        ax.set_xticks([]); ax.set_yticks([])
        t = f"slab {slab}"
        if per_slab_text:
            t += f"   {per_slab_text[slab]}"
        ax.set_title(t, fontsize=12)
    cb = fig.colorbar(im, ax=axes, shrink=0.7, pad=0.01)
    cb.set_label("log2( hits/event sim / data )   red = simulation has more", fontsize=12)
    fig.suptitle(title, fontsize=15)
    fig.savefig(path, dpi=110); plt.close(fig)


def local_residual(lr, ix, iy, half=2):
    """lr minus the median of its (2*half+1)^2 neighbourhood in the same slab: removes the smooth shower-shape and
    beam-position pattern and leaves what is local to a cell or a chip."""
    out = np.full(NCELL, np.nan)
    for slab in range(N_SLAB):
        g = slab_grid(lr, ix, iy, slab)
        sm = np.full_like(g, np.nan)
        for j in range(32):
            for i in range(32):
                w = g[max(j - half, 0):j + half + 1, max(i - half, 0):i + half + 1]
                if np.isfinite(g[j, i]) and np.isfinite(w).sum() >= 6:
                    sm[j, i] = np.nanmedian(w)
        sl = slice(slab * N_CHIP * N_CHAN, (slab + 1) * N_CHIP * N_CHAN)
        x, y = ix[sl], iy[sl]; k = x >= 0
        r = np.full(N_CHIP * N_CHAN, np.nan)
        r[k] = g[y[k], x[k]] - sm[y[k], x[k]]
        out[sl] = r
    return out


def group_variance(lr, err, groups):
    """Fraction of the scatter of lr (above its statistical floor) that the group means explain, one slab."""
    k = np.isfinite(lr)
    lr, err, groups = lr[k], err[k], groups[k]
    if len(lr) < 20:
        return np.nan, np.nan, np.nan
    w = 1 / err ** 2
    mu = np.average(lr, weights=w)
    tot = np.average((lr - mu) ** 2, weights=w)
    gm = {g: np.average(lr[groups == g], weights=w[groups == g]) for g in np.unique(groups)}
    within = np.average((lr - np.array([gm[g] for g in groups])) ** 2, weights=w)
    stat = np.average(err ** 2, weights=w)
    excess = max(tot - stat, 1e-12)
    return float(np.sqrt(tot)), float(np.sqrt(stat)), float(np.clip((tot - within) / excess, 0, 1.5))


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", required=True)
    ap.add_argument("--label", required=True)
    ap.add_argument("--data", required=True, help="data valtree(s), comma-separated")
    ap.add_argument("--sim", action="append", required=True, help="<name>=<valtree>")
    ap.add_argument("--trigger-table", help="Trigger_channel_highgain.txt of the set")
    ap.add_argument("--max-events", type=int, default=0, help="shower-like events per sample (0 = all)")
    ap.add_argument("--lim", type=float, default=1.0, help="colour range of the maps, |log2 ratio|")
    a = ap.parse_args(argv)
    os.makedirs(a.out, exist_ok=True)
    tag = a.label.replace(" ", "_").replace("/", "")

    d = load(a.data.split(","), a.max_events)
    sims = []
    for s in a.sim:
        name, path = s.split("=", 1)
        sims.append((name, load([path], a.max_events)))
    print(f"data: {d['nev']} shower-like events; " + ", ".join(f"{n}: {s['nev']}" for n, s in sims))
    ix, iy = pad_index(d["pos"])
    # cells the data never fired: take the position from the simulation
    miss = ix < 0
    if miss.any():
        ix2, iy2 = pad_index(sims[0][1]["pos"])
        ix[miss], iy[miss] = ix2[miss], iy2[miss]

    mu_cell = None
    if a.trigger_table:
        mu_cell, _ = read_trigger(a.trigger_table)

    chip_of = np.arange(NCELL) // N_CHAN
    # chip-blind control: the same 8x8 blocks shifted by 4 pads
    shifted = ((iy + 4) // 8) * 5 + (ix + 4) // 8 + 100 * (np.arange(NCELL) // (N_CHIP * N_CHAN))
    lines = [f"# {a.label}: data {d['nev']} shower-like events (" + a.data + ")"]
    coh = {}
    for si, (name, s) in enumerate(sims):
        lr, err, ok = ratio_cells(d, s)
        lrc = ratio_chips(d, s)
        lines.append(f"\n## {name}: {s['nev']} events; per slab: hits/event sim/data, cells drawn, "
                     f"rms log2 ratio, stat floor, chip share, shifted-block share; local residual: rms, chip share, shifted share")
        loc = local_residual(lr, ix, iy)
        per_slab = []
        coh[name] = []
        for slab in range(N_SLAB):
            sl = slice(slab * N_CHIP * N_CHAN, (slab + 1) * N_CHIP * N_CHAN)
            rs = (s["counts"][sl].sum() / s["nev"]) / max(d["counts"][sl].sum() / d["nev"], 1e-12)
            rms, stat, fchip = group_variance(lr[sl], err[sl], chip_of[sl])
            _, _, fshift = group_variance(lr[sl], err[sl], shifted[sl])
            lrms, _, lchip = group_variance(loc[sl], err[sl], chip_of[sl])
            _, _, lshift = group_variance(loc[sl], err[sl], shifted[sl])
            coh[name].append((rms, stat, fchip, fshift, lrms, lchip, lshift))
            per_slab.append(f"sim/data {rs:.3f}")
            lines.append(f"slab {slab:2d}  {rs:6.3f}  {int(ok[sl].sum()):5d}  {rms:6.3f}  {stat:6.3f}  "
                         f"{fchip:5.2f}  {fshift:5.2f}   {lrms:6.3f}  {lchip:5.2f}  {lshift:5.2f}")
        tot = s["counts"].sum() / s["nev"] / (d["counts"].sum() / d["nev"])
        lines.append(f"all slabs: hits/event sim/data {tot:.3f}")
        sname = name.replace(" ", "_")
        map_figure(lr, ix, iy, f"{a.label} - hits per cell, simulation ({name}) / data", f"{a.out}/cells_{tag}_{sname}.png",
                   a.lim, per_slab)
        map_figure(loc, ix, iy, f"{a.label} - hits per cell, simulation ({name}) / data, minus the median of the 5x5 "
                   "neighbourhood (local structure only)", f"{a.out}/cellsloc_{tag}_{sname}.png", 0.6, per_slab)
        map_figure(np.repeat(lrc, N_CHAN), ix, iy, f"{a.label} - hits per chip, simulation ({name}) / data",
                   f"{a.out}/chips_{tag}_{sname}.png", a.lim, per_slab)
        lines.append("chips (slab chip log2 ratio): " + " ".join(
            f"{c // N_CHIP}:{c % N_CHIP}:{v:+.2f}" for c, v in enumerate(lrc) if np.isfinite(v)))

    # coherence figure
    ncol = 3 if mu_cell is not None else 2
    fig, axes = plt.subplots(1, ncol, figsize=(8.5 * ncol, 6.2), squeeze=False, constrained_layout=True)
    ax, ax2 = axes[0, 0], axes[0, 1]
    slabs = np.arange(N_SLAB)
    for si, (name, _) in enumerate(sims):
        c = np.array(coh[name])
        ax.plot(slabs + 0.08 * si, c[:, 2], "o-", color=C_SIM[si], lw=2, ms=8, label=f"{name}: chips")
        ax.plot(slabs + 0.08 * si, c[:, 3], "s:", color=C_SIM[si], lw=1.6, ms=7, mfc="white",
                label=f"{name}: 8x8 blocks shifted 4 pads (chip-blind)")
        ax2.plot(slabs + 0.08 * si, c[:, 5], "o-", color=C_SIM[si], lw=2, ms=8, label=f"{name}: chips")
        ax2.plot(slabs + 0.08 * si, c[:, 6], "s:", color=C_SIM[si], lw=1.6, ms=7, mfc="white",
                 label=f"{name}: shifted 8x8 blocks (chip-blind)")
    ax.axhline(0, color="0.5", lw=1); ax.axhline(1, color="0.5", lw=1, ls="--")
    ax.set_xlabel("slab", fontsize=13); ax.set_ylabel("share of the cell scatter above statistics\nexplained by the group means", fontsize=13)
    ax.set_xticks(slabs); ax.grid(color=C_GRID, lw=0.6); ax.legend(fontsize=10, frameon=False)
    ax.set_title("Raw ratio: is the excess coherent within a chip?", fontsize=13)
    ax2.axhline(0, color="0.5", lw=1); ax2.axhline(1, color="0.5", lw=1, ls="--")
    ax2.set_xlabel("slab", fontsize=13); ax2.set_ylabel("share of the local-residual scatter above statistics\nexplained by the group means", fontsize=13)
    ax2.set_xticks(slabs); ax2.grid(color=C_GRID, lw=0.6); ax2.legend(fontsize=10, frameon=False)
    ax2.set_title("Local residual (smooth pattern removed)", fontsize=13)
    if mu_cell is not None:
        ax = axes[0, 2]
        chip_mu = np.array([np.nanmean(mu_cell[c * N_CHAN:(c + 1) * N_CHAN]) if np.isfinite(mu_cell[c * N_CHAN:(c + 1) * N_CHAN]).any()
                            else np.nan for c in range(N_SLAB * N_CHIP)])
        dmu = mu_cell - np.repeat(chip_mu, N_CHAN)
        bins = np.arange(-6, 6.01, 1.0); bc = 0.5 * (bins[1:] + bins[:-1])
        lines.append("\n## log2(sim/data) per cell vs mu_adc(cell) - mean mu_adc(chip), weighted mean per bin")
        for si, (name, s) in enumerate(sims):
            lr, err, ok = ratio_cells(d, s)
            # subtract each chip's own mean ratio: only the channel-level part is left
            lrc = ratio_chips(d, s)
            res = lr - np.repeat(lrc, N_CHAN)
            k = np.isfinite(res) & np.isfinite(dmu)
            w = 1 / err[k] ** 2
            ib = np.digitize(dmu[k], bins) - 1
            m = np.array([np.average(res[k][ib == b], weights=w[ib == b]) if (ib == b).sum() > 10 else np.nan for b in range(len(bc))])
            e = np.array([1 / np.sqrt(w[ib == b].sum()) if (ib == b).sum() > 10 else np.nan for b in range(len(bc))])
            ax.errorbar(bc + 0.1 * si, m, yerr=e, fmt="o-", color=C_SIM[si], lw=2, ms=8, capsize=2, label=name)
            sl = np.polyfit(dmu[k], res[k], 1, w=np.sqrt(w))
            lines.append(f"{name}: slope {sl[0]:+.4f} log2/ADC; " + " ".join(f"{x:+.1f}:{y:+.3f}" for x, y in zip(bc, m)))
        ax.axhline(0, color="0.5", lw=1)
        ax.set_xlabel("muon trigger threshold: cell - its chip mean [ADC]", fontsize=13)
        ax.set_ylabel("log2(sim/data) per cell - its chip's", fontsize=13)
        ax.grid(color=C_GRID, lw=0.6); ax.legend(fontsize=11, frameon=False)
        ax.set_title("Channel-level excess against the channel's own discriminator", fontsize=13)
    fig.suptitle(a.label, fontsize=15)
    fig.savefig(f"{a.out}/coherence_{tag}.png", dpi=110); plt.close(fig)

    # chip occupancy: hits in a chip vs the energy it collects, and the near-threshold hit spectrum
    fig, axes = plt.subplots(1, 3, figsize=(23, 6.4), constrained_layout=True)
    ebins = np.array([0.5, 1, 2, 3, 5, 8, 12, 20, 35, 60, 100, 200, 400, 1000.0])
    ec = np.sqrt(ebins[1:] * ebins[:-1])
    prof = {}
    for name, s, col, mk in [("data", d, C_DATA, "o")] + [(n, s, C_SIM[i], "s") for i, (n, s) in enumerate(sims)]:
        ib = np.digitize(s["chip_e"], ebins) - 1
        m = np.array([s["chip_n"][ib == b].mean() if (ib == b).sum() > 20 else np.nan for b in range(len(ec))])
        prof[name] = m
        axes[0].plot(ec, m, mk + "-", color=col, lw=2, ms=8, label=name)
        h, eb = np.histogram(s["hit_e"], bins=np.arange(0.5, 4.01, 0.1))
        axes[2].step(eb[:-1], h / s["nev"], where="post", color=col, lw=2, label=name)
    for i, (name, _) in enumerate(sims):
        axes[1].plot(ec, prof[name] / prof["data"], "s-", color=C_SIM[i], lw=2, ms=8, label=name)
    axes[1].axhline(1, color="0.5", lw=1, ls="--")
    for ax in axes[:2]:
        ax.set_xscale("log"); ax.set_xlabel("energy the chip collects in the event [MIP]", fontsize=13)
        ax.grid(color=C_GRID, lw=0.6); ax.legend(fontsize=11, frameon=False)
    axes[0].set_yscale("log"); axes[0].set_ylabel("hits in the chip (mean)", fontsize=13)
    axes[1].set_ylabel("hits in the chip, sim / data", fontsize=13)
    axes[2].set_xlabel("hit energy [MIP]", fontsize=13); axes[2].set_ylabel("hits per event per 0.1 MIP", fontsize=13)
    axes[2].grid(color=C_GRID, lw=0.6); axes[2].legend(fontsize=11, frameon=False)
    axes[0].set_title("Chip multiplicity vs chip energy", fontsize=13)
    axes[1].set_title("At the same chip energy, does the simulation fire more cells?", fontsize=13)
    axes[2].set_title("Hit energies near the threshold", fontsize=13)
    fig.suptitle(a.label, fontsize=15)
    fig.savefig(f"{a.out}/chipocc_{tag}.png", dpi=110); plt.close(fig)
    lines.append("\n## mean hits in a chip vs chip energy [MIP] (bin centres " + " ".join(f"{x:.3g}" for x in ec) + ")")
    for name, m in prof.items():
        lines.append(f"{name}: " + " ".join(f"{x:.2f}" for x in m))

    with open(f"{a.out}/trigger_maps_{tag}.txt", "w") as f:
        f.write("\n".join(lines) + "\n")
    print("\n".join(lines[:40]))


if __name__ == "__main__":
    main()

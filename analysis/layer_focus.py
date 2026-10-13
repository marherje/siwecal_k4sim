#!/usr/bin/env python
"""
Calibration effects behind the sim/data hit excess, per layer and for one layer in detail. Shower-like events of one
beam group from the valtrees (0.5 MIP hit cut in data and simulation), one or more digitisations of the same sample.

    python analysis/layer_focus.py --out <dir> --label "th210 S1 20 GeV" --slab 6 \
        --data <data valtree> --sim "per-slab=<valtree>" --sim "muon tables=<valtree>" \
        --calib /eos/experiment/drdcalo/siw-ecal/TB2026-06/Calib/MuonCalib --set th210

Figures (<tag> = label):
  energy_<tag>.png    sim/data hits against hit energy, all layers and the focus slab, for the cells whose muon
                      trigger threshold sits > 1.5 ADC above their chip's, the cells within 1 ADC of it, and the cells
                      whose local residual (minus the 5x5 median) is > +0.4 in log2: a threshold shift only removes
                      hits near the threshold, an inefficiency removes them at every energy
  layers_<tag>.png    per layer: hits/event sim/data in three hit-energy ranges, and the most probable hit energy
                      (peak of the 0.6-2.5 MIP spectrum) in data and simulation
  slab<N>_<tag>.png   the focus slab: data hits/event, sim/data for each simulation, local residual, muon trigger
                      threshold, its fit status, MIP MPV and pedestal width per cell
  spectrum<N>_<tag>.png  hit-energy spectra of the focus slab and its neighbours
  layer_focus_<tag>.txt  the numbers
"""
import argparse
import glob
import os

import numpy as np
import uproot
import awkward as ak
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.colors import TwoSlopeNorm  # noqa: E402

from trigger_maps import (N_SLAB, N_CHIP, N_CHAN, NCELL, C_DATA, C_SIM, C_GRID, MIN_HITS, showerlike, pad_index,
                          slab_grid, chip_outlines, local_residual, read_trigger)

EBINS = np.array([0.5, 0.8, 1.2, 2.0, 4.0, 10.0, 1e4])
ELAB = ["0.5-0.8", "0.8-1.2", "1.2-2", "2-4", "4-10", ">10"]
SPEC = np.arange(0.5, 6.01, 0.1)


def load(path, nmax=0):
    """Per-cell hit counts in EBINS, per-slab spectra, cell positions."""
    cnt = np.zeros((NCELL, len(EBINS) - 1))
    spec = np.zeros((N_SLAB, len(SPEC) - 1))
    xs, ys, nxy = np.zeros(NCELL), np.zeros(NCELL), np.zeros(NCELL)
    nev = 0
    for a in uproot.iterate(f"{path}:ecal", ["nhit_chan", "hit_slab", "hit_chip", "hit_chan", "hit_x", "hit_y",
                                             "hit_energy"], step_size="300 MB", library="ak"):
        a = a[showerlike(ak.to_numpy(a["nhit_chan"]).astype(float))]
        if nmax and nev + len(a) > nmax:
            a = a[: nmax - nev]
        sl = ak.to_numpy(ak.flatten(a["hit_slab"])).astype(int)
        cell = (sl * N_CHIP + ak.to_numpy(ak.flatten(a["hit_chip"])).astype(int)) * N_CHAN \
            + ak.to_numpy(ak.flatten(a["hit_chan"])).astype(int)
        e = ak.to_numpy(ak.flatten(a["hit_energy"])).astype(float)
        ib = np.clip(np.digitize(e, EBINS) - 1, 0, len(EBINS) - 2)
        np.add.at(cnt, (cell, ib), 1)
        k = (e >= SPEC[0]) & (e < SPEC[-1])
        np.add.at(spec, (sl[k], ((e[k] - SPEC[0]) / 0.1).astype(int)), 1)
        xs += np.bincount(cell, ak.to_numpy(ak.flatten(a["hit_x"])), minlength=NCELL)
        ys += np.bincount(cell, ak.to_numpy(ak.flatten(a["hit_y"])), minlength=NCELL)
        nxy += np.bincount(cell, minlength=NCELL)
        nev += len(a)
        if nmax and nev >= nmax:
            break
    with np.errstate(invalid="ignore", divide="ignore"):
        pos = np.stack([xs / nxy, ys / nxy], axis=1)
    return {"cnt": cnt, "spec": spec, "nev": nev, "pos": pos}


def read_cells(path, col, bad=lambda f: False):
    """One column per (layer, chip, channel) of a calibration table, NaN where bad."""
    v = np.full(NCELL, np.nan)
    for line in open(path):
        if line.startswith("#"):
            continue
        f = line.split()
        if len(f) <= col:
            continue
        c = (int(f[0]) * N_CHIP + int(f[1])) * N_CHAN + int(f[2])
        if not bad(f):
            v[c] = float(f[col])
    return v


def ratio(d, s, cells, ebin=None):
    """hits/event sim/data summed over cells (and an energy bin), with its error."""
    dc = d["cnt"][cells] if ebin is None else d["cnt"][cells][:, ebin]
    sc = s["cnt"][cells] if ebin is None else s["cnt"][cells][:, ebin]
    nd, ns = dc.sum(), sc.sum()
    if nd < MIN_HITS or ns < MIN_HITS:
        return np.nan, np.nan
    r = (ns / s["nev"]) / (nd / d["nev"])
    return r, r * np.sqrt(1 / nd + 1 / ns)


def peak(h):
    """Most probable value of a 0.1-MIP spectrum in 0.6-2.5 MIP: parabola through the maximum bin."""
    c = 0.5 * (SPEC[1:] + SPEC[:-1])
    w = (c > 0.6) & (c < 2.5)
    i = np.flatnonzero(w)[np.argmax(h[w])]
    if 0 < i < len(h) - 1:
        y0, y1, y2 = h[i - 1], h[i], h[i + 1]
        den = y0 - 2 * y1 + y2
        return c[i] + (0.05 * (y0 - y2) / den if den != 0 else 0)
    return c[i]


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", required=True)
    ap.add_argument("--label", required=True)
    ap.add_argument("--slab", type=int, default=6)
    ap.add_argument("--data", required=True)
    ap.add_argument("--sim", action="append", required=True)
    ap.add_argument("--calib", required=True, help="calibration dir with trigger/, mips/, pedestals/")
    ap.add_argument("--set", default="th210")
    ap.add_argument("--max-events", type=int, default=0)
    a = ap.parse_args(argv)
    os.makedirs(a.out, exist_ok=True)
    tag = a.label.replace(" ", "_")
    S = a.slab

    d = load(a.data, a.max_events)
    sims = [(n, load(p, a.max_events)) for n, p in (s.split("=", 1) for s in a.sim)]
    ix, iy = pad_index(d["pos"])
    miss = ix < 0
    ix2, iy2 = pad_index(sims[0][1]["pos"])
    ix[miss], iy[miss] = ix2[miss], iy2[miss]

    mu, st = read_trigger(f"{a.calib}/trigger/{a.set}/Trigger_channel_highgain.txt")
    chip_mu = np.array([np.nanmean(mu[c * N_CHAN:(c + 1) * N_CHAN]) if np.isfinite(mu[c * N_CHAN:(c + 1) * N_CHAN]).any()
                        else np.nan for c in range(N_SLAB * N_CHIP)])
    dmu = mu - np.repeat(chip_mu, N_CHAN)
    mip_file = glob.glob(f"{a.calib}/mips/{a.set}/MIP_*_highgain.txt")[0]
    mpv = read_cells(mip_file, 3, bad=lambda f: float(f[3]) <= 0)
    mpv_own = read_cells(mip_file, 4, bad=lambda f: float(f[4]) <= 0)   # empv > 0: the channel's own fit
    ped_file = glob.glob(f"{a.calib}/pedestals/{a.set}/Pedestal_*_highgain.txt")[0]
    # median over the SCAs of the pedestal width (columns 5, 8, 11, ...; error < 0 = no fit)
    pw = np.full(NCELL, np.nan)
    for line in open(ped_file):
        if line.startswith("#"):
            continue
        f = line.split()
        c = (int(f[0]) * N_CHIP + int(f[1])) * N_CHAN + int(f[2])
        w = [float(f[3 + 3 * k + 2]) for k in range((len(f) - 3) // 3) if float(f[3 + 3 * k + 1]) >= 0]
        w = [x for x in w if x > 0]
        if w:
            pw[c] = np.median(w)

    # local residual of the first simulation, total hits
    dt = {"counts": d["cnt"].sum(1), "nev": d["nev"]}
    lines = [f"# {a.label}: data {d['nev']} shower-like events; " + ", ".join(f"{n} {s['nev']}" for n, s in sims)]
    lr_all, loc_all = {}, {}
    for n, s in sims:
        st_ = {"counts": s["cnt"].sum(1), "nev": s["nev"]}
        ok = (dt["counts"] >= MIN_HITS) & (st_["counts"] >= MIN_HITS)
        with np.errstate(divide="ignore", invalid="ignore"):
            lr = np.log2((st_["counts"] / s["nev"]) / (dt["counts"] / d["nev"]))
        lr[~ok] = np.nan
        lr_all[n] = lr
        loc_all[n] = local_residual(lr, ix, iy)
    first = sims[0][0]
    slab_of = np.arange(NCELL) // (N_CHIP * N_CHAN)
    classes = [("muon threshold > chip + 1.5 ADC", np.isfinite(dmu) & (dmu > 1.5)),
               ("muon threshold within 1 ADC of chip", np.isfinite(dmu) & (np.abs(dmu) < 1.0)),
               (f"local residual > +0.4 ({first})", np.nan_to_num(loc_all[first], nan=-9) > 0.4)]

    # 1) energy-binned diagnosis
    fig, axes = plt.subplots(1, 2, figsize=(17, 6.4), constrained_layout=True, sharey=True)
    lines.append("\n## sim/data hits per hit-energy bin [MIP] " + " ".join(ELAB))
    for ax, (scope, sel) in zip(axes, [("all layers", np.ones(NCELL, bool)), (f"slab {S}", slab_of == S)]):
        for ci, (cname, cm) in enumerate(classes):
            for si, (n, s) in enumerate(sims):
                r = np.array([ratio(d, s, cm & sel, b) for b in range(len(ELAB))])
                ls = ["-", "--", ":"][ci]
                mk = ["o", "s", "^"][ci]
                ax.errorbar(np.arange(len(ELAB)) + 0.06 * (si + 2 * ci) - 0.2, r[:, 0], yerr=r[:, 1], fmt=mk + ls,
                            color=C_SIM[si], lw=1.8, ms=7, capsize=2, mfc="white" if ci == 1 else C_SIM[si],
                            label=f"{cname} - {n}")
                lines.append(f"{scope} | {cname} | {n} ({int((cm & sel).sum())} cells): "
                             + " ".join(f"{x:.3f}" for x in r[:, 0]))
        ax.axhline(1, color="0.5", lw=1, ls="--")
        ax.set_xticks(range(len(ELAB))); ax.set_xticklabels(ELAB)
        ax.set_xlabel("hit energy [MIP]", fontsize=13); ax.grid(color=C_GRID, lw=0.6)
        ax.set_title(f"{scope}: sim/data hits per event by hit energy", fontsize=13)
    axes[0].set_ylabel("hits per event, simulation / data", fontsize=13)
    axes[1].legend(fontsize=9, frameon=False, loc="upper right")
    fig.suptitle(a.label, fontsize=15)
    fig.savefig(f"{a.out}/energy_{tag}.png", dpi=110); plt.close(fig)

    # 2) per layer
    fig, axes = plt.subplots(1, 2, figsize=(17, 6.2), constrained_layout=True)
    rng = [("< 1.2 MIP", [0, 1]), ("1.2-4 MIP", [2, 3]), ("> 4 MIP", [4, 5])]
    lines.append("\n## per layer: sim/data hits/event (<1.2, 1.2-4, >4 MIP, all); peak data / sims")
    pk = {"data": [peak(d["spec"][L]) for L in range(N_SLAB)]}
    for n, s in sims:
        pk[n] = [peak(s["spec"][L]) for L in range(N_SLAB)]
    for si, (n, s) in enumerate(sims):
        for ri, (rl, bins) in enumerate(rng):
            r = []
            for L in range(N_SLAB):
                cm = slab_of == L
                nd = d["cnt"][cm][:, bins].sum(); ns = s["cnt"][cm][:, bins].sum()
                r.append((ns / s["nev"]) / max(nd / d["nev"], 1e-12))
            axes[0].plot(np.arange(N_SLAB) + 0.07 * si, r, ["o-", "s--", "^:"][ri], color=C_SIM[si], lw=1.8, ms=7,
                         label=f"{n}: {rl}")
    for L in range(N_SLAB):
        cm = slab_of == L
        vals = []
        for n, s in sims:
            vals += [(s["cnt"][cm][:, b].sum() / s["nev"]) / max(d["cnt"][cm][:, b].sum() / d["nev"], 1e-12)
                     for b in ([0, 1], [2, 3], [4, 5])]
            vals.append((s["cnt"][cm].sum() / s["nev"]) / max(d["cnt"][cm].sum() / d["nev"], 1e-12))
        lines.append(f"layer {L:2d}: " + " ".join(f"{v:.3f}" for v in vals) + "   peak " +
                     " ".join(f"{pk[k][L]:.3f}" for k in pk))
    axes[0].axhline(1, color="0.5", lw=1, ls="--")
    axes[0].axvspan(S - 0.4, S + 0.4, color=C_GRID, alpha=0.6, zorder=0)
    axes[0].set_xlabel("layer", fontsize=13); axes[0].set_ylabel("hits per event, simulation / data", fontsize=13)
    axes[0].set_xticks(range(N_SLAB)); axes[0].grid(color=C_GRID, lw=0.6); axes[0].legend(fontsize=9, frameon=False, ncol=2)
    axes[0].set_title("Excess per layer, split by hit energy", fontsize=13)
    axes[1].plot(range(N_SLAB), pk["data"], "o-", color=C_DATA, lw=2, ms=8, label="data")
    for si, (n, _) in enumerate(sims):
        axes[1].plot(np.arange(N_SLAB) + 0.07 * (si + 1), pk[n], "s--", color=C_SIM[si], lw=1.8, ms=7, label=n)
    axes[1].axvspan(S - 0.4, S + 0.4, color=C_GRID, alpha=0.6, zorder=0)
    axes[1].set_xlabel("layer", fontsize=13); axes[1].set_ylabel("most probable hit energy [MIP]", fontsize=13)
    axes[1].set_xticks(range(N_SLAB)); axes[1].grid(color=C_GRID, lw=0.6); axes[1].legend(fontsize=11, frameon=False)
    axes[1].set_title("Energy scale per layer: peak of the hit spectrum", fontsize=13)
    fig.suptitle(a.label, fontsize=15)
    fig.savefig(f"{a.out}/layers_{tag}.png", dpi=110); plt.close(fig)

    # 3) focus slab maps
    panels = [("data: hits per event", np.log10(np.where(dt["counts"] > 0, dt["counts"] / d["nev"], np.nan)), "viridis",
               None, "log10 hits/event")]
    for n, _ in sims:
        panels.append((f"{n}: log2 sim/data", lr_all[n], "RdBu_r", (-1, 1), "log2"))
    panels.append((f"{sims[-1][0]}: local residual", loc_all[sims[-1][0]], "RdBu_r", (-0.6, 0.6), "log2"))
    panels.append(("muon trigger threshold [ADC]", mu, "magma", None, "ADC"))
    panels.append(("muon threshold - chip mean [ADC]", dmu, "PuOr_r", (-3, 3), "ADC"))
    panels.append(("MIP MPV [ADC] (own fits only)", np.where(np.isfinite(mpv_own), mpv, np.nan), "cividis", None, "ADC"))
    panels.append(("pedestal width [ADC] (median over SCAs)", pw, "cividis", None, "ADC"))
    nc = 4
    nr = int(np.ceil(len(panels) / nc))
    fig, axes = plt.subplots(nr, nc, figsize=(5.4 * nc, 5.0 * nr), constrained_layout=True)
    for ax, (t, v, cm, lim, cl) in zip(axes.ravel(), panels):
        g = slab_grid(v, ix, iy, S)
        cmap = plt.get_cmap(cm).copy(); cmap.set_bad("#bdbdbd")
        if lim:
            im = ax.imshow(g, origin="lower", cmap=cmap, norm=TwoSlopeNorm(vcenter=0, vmin=lim[0], vmax=lim[1]),
                           interpolation="nearest")
        else:
            f = g[np.isfinite(g)]
            im = ax.imshow(g, origin="lower", cmap=cmap, interpolation="nearest",
                           vmin=np.percentile(f, 2) if f.size else None, vmax=np.percentile(f, 98) if f.size else None)
        chip_outlines(ax, ix, iy, S)
        ax.set_xticks([]); ax.set_yticks([]); ax.set_title(t, fontsize=12)
        fig.colorbar(im, ax=ax, shrink=0.8, label=cl)
    for ax in axes.ravel()[len(panels):]:
        ax.axis("off")
    fig.suptitle(f"{a.label} - slab {S}", fontsize=15)
    fig.savefig(f"{a.out}/slab{S}_{tag}.png", dpi=105); plt.close(fig)

    # slab numbers: calibration quality per slab
    lines.append("\n## per slab: muon trigger own fits / chip / masked, MIP own fits, median MPV, median ped width, "
                 "median |dmu|")
    for L in range(N_SLAB):
        cm = slab_of == L
        lines.append(f"slab {L:2d}: trig {int((st[cm] == 0).sum())}/{int((st[cm] == 1).sum())}/{int((st[cm] == 3).sum())}"
                     f"  mip own {int(np.isfinite(mpv_own[cm]).sum())}  MPV {np.nanmedian(mpv[cm]):.2f}"
                     f"  ped w {np.nanmedian(pw[cm]):.2f}  chip mu {np.nanmedian(chip_mu[L * N_CHIP:(L + 1) * N_CHIP]):.2f}"
                     f"  |dmu| {np.nanmedian(np.abs(dmu[cm])):.2f}")
    # hot cells of the focus slab
    hot = np.flatnonzero((slab_of == S) & (np.nan_to_num(loc_all[first], nan=-9) > 0.4))
    lines.append(f"\n## slab {S}: cells with local residual > 0.4 ({first}): chip chan dmu mu MPV(own) pedw "
                 "data_hits/ev " + " ".join(f"loc[{n}]" for n, _ in sims))
    for c in hot:
        lines.append(f"  {c // N_CHAN % N_CHIP:2d} {c % N_CHAN:2d}  {dmu[c]:+.2f}  {mu[c]:.2f}  {mpv_own[c] if np.isfinite(mpv_own[c]) else np.nan:.3f}"
                     f"  {pw[c]:.2f}  {dt['counts'][c] / d['nev']:.4f}  " + " ".join(f"{loc_all[n][c]:+.2f}" for n, _ in sims))

    # 4) spectra of the focus slab and neighbours
    nb = [L for L in (S - 1, S, S + 1) if 0 <= L < N_SLAB]
    fig, axes = plt.subplots(1, len(nb), figsize=(7.2 * len(nb), 5.6), constrained_layout=True, sharey=True)
    c = 0.5 * (SPEC[1:] + SPEC[:-1])
    for ax, L in zip(np.atleast_1d(axes), nb):
        ax.step(c, d["spec"][L] / d["nev"], where="mid", color=C_DATA, lw=2.2, label=f"data (peak {pk['data'][L]:.2f})")
        for si, (n, s) in enumerate(sims):
            ax.step(c, s["spec"][L] / s["nev"], where="mid", color=C_SIM[si], lw=1.8, ls=["-", "--"][si % 2],
                    label=f"{n} (peak {pk[n][L]:.2f})")
        ax.set_title(f"slab {L}" + ("  (focus)" if L == S else ""), fontsize=13)
        ax.set_xlabel("hit energy [MIP]", fontsize=13); ax.grid(color=C_GRID, lw=0.6); ax.legend(fontsize=10, frameon=False)
    np.atleast_1d(axes)[0].set_ylabel("hits per event per 0.1 MIP", fontsize=13)
    fig.suptitle(a.label, fontsize=15)
    fig.savefig(f"{a.out}/spectrum{S}_{tag}.png", dpi=110); plt.close(fig)

    with open(f"{a.out}/layer_focus_{tag}.txt", "w") as f:
        f.write("\n".join(lines) + "\n")
    print("\n".join(lines))


if __name__ == "__main__":
    main()

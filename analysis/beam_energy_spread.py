#!/usr/bin/env python3
"""Beam energy spread scan: which fractional Gaussian spread sigma_E of the beam momentum makes the digitised
simulation reproduce the calibrated-energy distribution of the data, per energy, on the th230 energy scans at E1 and
E2 (the only positions where the resolution is studied, see docs/plan_border_cut_final_v8.md).

Input: the per-event caches of final_panels.py (sum_w of data, nodigi and digi after shower-like + border cut), E1
with a 45 mm border cut and E2 with 40 mm, calibrated with the th230 MIP -> GeV line of each sample. The simulation
(final_v9) was generated with sigma_E0 = 2 % and keeps no per-event truth energy, so other spreads are emulated:

  sigma_E >= sigma_E0  smearing: every simulated event gets E -> E * (1 + d),
                       d ~ N(0, sqrt(sigma_E^2 - sigma_E0^2)) (exact for a linear response; common random numbers
                       across the grid)
  sigma_E <  sigma_E0  width scaling: E -> mu + k (E - mu) around the Gaussian-core mean, with
                       k = sqrt(r0^2 - sigma_E0^2 + sigma_E^2) / r0 and r0 the core sigma/mu of the 2 % sample
                       (the beam term adds in quadrature to the calorimeter one). Approximation: it also scales the
                       tails; validated against the smearing above 2 % (--closure) and, once it exists, against a
                       sample generated at 1 % (--sim-1pc).

Two figures of merit per point and sigma_E:
  chi2_hist  binned chi2 of the unit-area energy distributions in the Gaussian-core window of the data
             [mu - 1.5 sigma, mu + 2 sigma] (each sample divided by its own core mean, so only the shape is compared)
  chi2_res   (r_data - r_sim)^2 / (err_data^2 + err_sim^2) on the core sigma/mu, err = r / sqrt(2 N)
The best sigma_E per point is the minimum of a parabola through the chi2_hist scan; its error is chi2_min + 1.

    python3 analysis/beam_energy_spread.py <out>  [--grid 1.0:3.0:0.2] [--sim-1pc <dir of 1 % valtrees>]
"""
import argparse
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from analysis import plot_style as ps  # noqa: E402
from analysis.compare_event_level import gauss_core  # noqa: E402

B = "/eos/experiment/drdcalo/siw-ecal/TB2026-06/Simulation/tests/border_cut"
# (scan, border cut, final_panels outdir with the cache and the MIP -> GeV line)
SCANS = {"E1": (45, f"{B}/figs_v9_p9/b45"), "E2": (40, f"{B}/figs_v9_setcal/b40")}
ENERGIES = [7.5, 10, 20, 34, 52, 74]
# Points used in the fits: only runs taken after the XCSH collimator was closed (+-30 -> +-20 mm, around run 20); the
# earlier ones (E1 20-74 GeV, runs 12-20) have a dirtier beam, and E1 52/74 GeV no clean selection at all.
# Filled by main() from the runs column of beam_points.txt and --min-run.
NOFIT = set()
POINTS = os.path.join(os.path.dirname(os.path.abspath(__file__)), "beam_points.txt")


def point_runs(scan, E):
    for line in open(POINTS):
        f = line.split("#")[0].split()
        if len(f) >= 9 and f[0] == f"th230_{scan}_e{E:g}":
            return [int(r) for r in f[8].split(",")]
    raise SystemExit(f"{POINTS}: no th230_{scan}_e{E:g}")
SIGMA_E0 = 0.02      # final_v9 generation
H2_MAX = 0.02        # quoted H2 maximum
NBINS = 30


def th230_line(path, kind):
    for line in open(path):
        f = line.split("#")[0].split()
        if len(f) == 3 and f[0] == f"th230:{kind}":
            return float(f[1]), float(f[2])
    raise SystemExit(f"{path}: no th230:{kind} line")


def load(scan, E, kinds=("data", "digi")):
    cut, d = SCANS[scan]
    z = np.load(f"{d}/cache_final_v9_b{cut}/th230_{scan}_e{E:g}.npz")
    out = {}
    for k in kinds:
        a, b = th230_line(f"{d}/mip_gev_line.txt", k)
        out[k] = (z[f"{k}__sum_w"] - b) / a
    return out


def core(x):
    mu, sg = gauss_core(x)
    return mu, sg / mu


def emulate(x, sigma_e, rng_d, r0=None, mu0=None):
    """The 2 % sample turned into one with fractional beam spread sigma_e (see the module docstring)."""
    if mu0 is None:
        mu0, r0 = core(x)
    if sigma_e >= SIGMA_E0:
        return x * (1 + np.sqrt(sigma_e ** 2 - SIGMA_E0 ** 2) * rng_d)
    k = np.sqrt(max(r0 ** 2 - SIGMA_E0 ** 2 + sigma_e ** 2, 1e-12)) / r0
    return mu0 + k * (x - mu0)


def hist_chi2(data, sim):
    """Shape chi2 in the data core window, both samples divided by their own core mean and normalised to unit area
    in the window."""
    md, rd = core(data)
    ms, _ = core(sim)
    u, v = data / md, sim / ms
    edges = np.linspace(1 - 1.5 * rd, 1 + 2 * rd, NBINS + 1)
    cd, _ = np.histogram(u, edges)
    cs, _ = np.histogram(v, edges)
    nd, ns = cd.sum(), cs.sum()
    pd, ps_ = cd / nd, cs / ns
    var = cd / nd ** 2 + cs / ns ** 2
    ok = var > 0
    return float((((pd - ps_) ** 2)[ok] / var[ok]).sum()), int(ok.sum()) - 1


def parabola_min(xs, ys):
    """Minimum of a parabola through the 3 lowest neighbouring points (x in %), and the chi2_min + 1 half-width."""
    i = int(np.argmin(ys))
    i = min(max(i, 1), len(xs) - 2)
    c = np.polyfit(xs[i - 1:i + 2], ys[i - 1:i + 2], 2)
    if c[0] <= 0:
        return xs[int(np.argmin(ys))], np.nan, float(np.min(ys)), False
    xm = -c[1] / (2 * c[0])
    ym = np.polyval(c, xm)
    inside = xs[0] <= xm <= xs[-1]
    return float(xm), float(1 / np.sqrt(c[0])), float(ym), inside


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("out")
    ap.add_argument("--grid", default="1.0:3.0:0.2", help="sigma_E grid in %%, lo:hi:step")
    ap.add_argument("--seed", type=int, default=12345)
    ap.add_argument("--min-run", type=int, default=21,
                    help="fit only points whose runs are all >= MIN_RUN (collimator closed after run 20)")
    args = ap.parse_args(argv)
    RUNS = {(sc, E): point_runs(sc, E) for sc in SCANS for E in ENERGIES}
    NOFIT.update(k for k, r in RUNS.items() if min(r) < args.min_run)
    print("left out of the fits (runs < %d): %s" % (args.min_run, sorted(NOFIT)))
    os.makedirs(args.out, exist_ok=True)
    lo, hi, st = (float(v) for v in args.grid.split(":"))
    grid = np.round(np.arange(lo, hi + st / 2, st), 3)          # %
    plt.rcParams.update(ps.RC)

    R = {}   # (scan, E) -> dict
    for scan in SCANS:
        for E in ENERGIES:
            s = load(scan, E)
            data, sim = s["data"], s["digi"]
            rng = np.random.default_rng(args.seed)
            dz = rng.standard_normal(len(sim))
            mu0, r0 = core(sim)
            md, rd = core(data)
            rows = []
            for g in grid:
                x = emulate(sim, g / 100, dz, r0, mu0)
                ch, ndf = hist_chi2(data, x)
                ms, rs = core(x)
                err = np.hypot(rd / np.sqrt(2 * len(data)), rs / np.sqrt(2 * len(x)))
                rows.append((g, ch, ndf, rs, ((rd - rs) / err) ** 2))
            rows = np.array(rows)
            # closure of the width scaling: above 2 %, scaling vs smearing
            clos = []
            for g in grid[grid > SIGMA_E0 * 100]:
                k = np.sqrt(r0 ** 2 - SIGMA_E0 ** 2 + (g / 100) ** 2) / r0
                xs_ = mu0 + k * (sim - mu0)
                clos.append((g, core(emulate(sim, g / 100, dz, r0, mu0))[1], core(xs_)[1]))
            bh = parabola_min(grid, rows[:, 1])
            br = parabola_min(grid, rows[:, 4])
            R[(scan, E)] = dict(rows=rows, rd=rd, r0=r0, Nd=len(data), Ns=len(sim), best_h=bh, best_r=br,
                                closure=np.array(clos), data=data, sim=sim, dz=dz, mu0=mu0)
            print(f"{scan} {E:5g} GeV  N {len(data)}/{len(sim)}  r_data {100*rd:.2f}%  r_sim(2%) {100*r0:.2f}%  "
                  f"best hist {bh[0]:.2f}+-{bh[1]:.2f}%  best res {br[0]:.2f}+-{br[1]:.2f}%")

    # ---------- text summary
    lines = ["# beam energy spread scan (final_v9 digi emulated, th230 E1 cut 45 mm / E2 cut 40 mm, th230 MIP->GeV line)",
             "# grid sigma_E [%]: " + " ".join(f"{g:.1f}" for g in grid),
             "# scan E N_data N_sim r_data r_sim(2%) | best_hist err chi2min/ndf | best_res err | chi2_hist/ndf per grid"]
    for (scan, E), r in R.items():
        bh, br = r["best_h"], r["best_r"]
        ndf = int(r["rows"][0, 2])
        flag = (" nofit" if (scan, E) in NOFIT else "") + " runs " + ",".join(map(str, RUNS[(scan, E)]))
        lines.append(f"{scan} {E:g} {r['Nd']} {r['Ns']} {100*r['rd']:.2f}% {100*r['r0']:.2f}% | "
                     f"{bh[0]:.2f} {bh[1]:.2f} {bh[2]:.1f}/{ndf} | {br[0]:.2f} {br[1]:.2f} | "
                     + " ".join(f"{v/ndf:.2f}" for v in r["rows"][:, 1]) + flag)
    # combined (clean points)
    tot_h = sum(r["rows"][:, 1] for k, r in R.items() if k not in NOFIT)
    tot_r = sum(r["rows"][:, 4] for k, r in R.items() if k not in NOFIT)
    ndf_h = sum(int(r["rows"][0, 2]) for k, r in R.items() if k not in NOFIT) - 1
    nclean = sum(1 for k in R if k not in NOFIT)
    ch = parabola_min(grid, tot_h)
    cr = parabola_min(grid, tot_r)
    in_range = (grid >= 1.0 - 1e-9) & (grid <= 2.0 + 1e-9)
    gi = grid[in_range]
    ch_in = gi[int(np.argmin(tot_h[in_range]))]
    lines += ["", f"combined over {nclean} clean points: chi2_hist min at {ch[0]:.2f} +- {ch[1]:.2f} % "
              f"(chi2/ndf {ch[2]:.0f}/{ndf_h}), chi2_res min at {cr[0]:.2f} +- {cr[1]:.2f} % "
              f"(chi2/ndf {cr[2]:.1f}/{nclean - 1}); best grid point in 1.0-2.0 %: {ch_in:.1f} %",
              "combined chi2_hist per grid: " + " ".join(f"{g:.1f}:{v:.0f}" for g, v in zip(grid, tot_h)),
              "combined chi2_res  per grid: " + " ".join(f"{g:.1f}:{v:.1f}" for g, v in zip(grid, tot_r)),
              "", "# closure (sigma_E > 2 %): core sigma/mu with smearing vs width scaling"]
    for (scan, E), r in R.items():
        lines.append(f"{scan} {E:g} " + " ".join(f"{g:.1f}:{100*a:.2f}/{100*b:.2f}" for g, a, b in r["closure"]))
    open(f"{args.out}/beam_energy_spread.txt", "w").write("\n".join(lines) + "\n")
    print("\n".join(lines[-(len(R) + 6):]))

    # ---------- figure 1: chi2_hist / ndf vs sigma_E, one panel per point
    fig, axs = plt.subplots(2, 6, figsize=(20, 7), sharex=True)
    for i, scan in enumerate(SCANS):
        for j, E in enumerate(ENERGIES):
            ax, r = axs[i, j], R[(scan, E)]
            ndf = r["rows"][0, 2]
            ax.axvspan(grid[0] - 0.1, 2.0, color="0.92", zorder=0)
            ax.plot(grid, r["rows"][:, 1] / ndf, "o-", color=ps.DIGI, ms=4, label="χ² shape / ndf")
            ax.plot(grid, r["rows"][:, 4] + r["rows"][:, 1].min() / ndf, "s--", color=ps.NODIGI, ms=3.5,
                    label="χ² σ/μ (+ offset)")
            bh = r["best_h"]
            if np.isfinite(bh[1]):
                ax.axvline(bh[0], color=ps.DIGI, lw=1)
                ax.axvspan(bh[0] - bh[1], bh[0] + bh[1], color=ps.DIGI, alpha=0.12)
            ax.axvline(100 * H2_MAX, color="k", ls=":", lw=1)
            ttl = f"{scan}  {E:g} GeV   best {bh[0]:.2f}" + (f" ± {bh[1]:.2f} %" if np.isfinite(bh[1]) else " %")
            ax.set_title(ttl + ("  (no fit)" if (scan, E) in NOFIT else ""), fontsize=9)
            if j == 0:
                ax.set_ylabel("χ² / ndf")
            if i == 1:
                ax.set_xlabel("beam σ_E [%]")
            if i == 0 and j == 0:
                ax.legend(fontsize=7)
    fig.suptitle("Beam energy spread scan: digitised sim (final_v9, emulated σ_E) vs data, calibrated energy, "
                 "th230 E1 (45 mm) / E2 (40 mm). Grey: 1–2 % (H2 range)", fontsize=10)
    fig.tight_layout()
    fig.savefig(f"{args.out}/chi2_per_point.png", dpi=130)
    plt.close(fig)

    # ---------- figure 2: best sigma_E vs energy + combined chi2
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(12, 4.6))
    for scan, mk, off in (("E1", "s", -0.6), ("E2", "o", 0.6)):
        for E in ENERGIES:
            r = R[(scan, E)]
            bh = r["best_h"]
            nf = (scan, E) in NOFIT
            a1.errorbar(E + off, bh[0], yerr=bh[1] if np.isfinite(bh[1]) else 0, fmt=mk, ms=6,
                        mfc="none" if nf else ps.DIGI, color=ps.DIGI, capsize=3,
                        label=f"{scan} (shape χ²)" if E == 7.5 else None)
            br = r["best_r"]
            a1.errorbar(E + off + 0.3 * np.sign(off), br[0], yerr=br[1] if np.isfinite(br[1]) else 0, fmt=mk,
                        ms=4.5, mfc="none" if nf else ps.NODIGI, color=ps.NODIGI, capsize=2,
                        label=f"{scan} (σ/μ χ²)" if E == 7.5 else None)
    a1.axhspan(1.0, 2.0, color="0.92", zorder=0, label="1–2 % (H2)")
    a1.axhline(ch[0], color="k", lw=1, label=f"combined: {ch[0]:.2f} %")
    a1.set_xlabel("beam energy [GeV]")
    a1.set_ylabel("best beam σ_E [%]")
    a1.set_ylim(0, 4)
    a1.legend(fontsize=7, ncol=2)
    a1.set_title("Best-fit spread per point (open: runs ≤ 20, not fitted)", fontsize=10)
    a2.plot(grid, tot_h - tot_h.min(), "o-", color=ps.DIGI, label=f"Σ shape χ² (ndf {ndf_h})")
    a2.plot(grid, tot_r - tot_r.min(), "s--", color=ps.NODIGI, label=f"Σ σ/μ χ² (ndf {nclean - 1})")
    a2.axvspan(grid[0] - 0.1, 2.0, color="0.92", zorder=0)
    a2.axvline(100 * H2_MAX, color="k", ls=":", lw=1)
    a2.set_yscale("symlog", linthresh=10)
    a2.set_ylim(bottom=0)
    a2.set_xlabel("beam σ_E [%]")
    a2.set_ylabel("Δχ² (combined, clean points)")
    a2.legend(fontsize=8)
    a2.set_title(f"Combined: shape {ch[0]:.2f} ± {ch[1]:.2f} %, σ/μ {cr[0]:.2f} ± {cr[1]:.2f} %", fontsize=10)
    fig.tight_layout()
    fig.savefig(f"{args.out}/best_sigmaE.png", dpi=130)
    plt.close(fig)

    # ---------- figure 3: sigma/mu vs E, data against the emulated sims
    fig, axs = plt.subplots(1, 2, figsize=(12, 4.6), sharey=True)
    show = [1.0, 1.4, 2.0, 2.6]
    cols = [st_['color'] for st_ in ps.reds(len(show))]
    for ax, scan in zip(axs, SCANS):
        Es = np.array(ENERGIES)
        rd = np.array([R[(scan, E)]["rd"] for E in ENERGIES])
        ed = rd / np.sqrt(2 * np.array([R[(scan, E)]["Nd"] for E in ENERGIES]))
        for g, c in zip(show, cols):
            gi_ = int(np.argmin(abs(grid - g)))
            rs = np.array([R[(scan, E)]["rows"][gi_, 3] for E in ENERGIES])
            ax.plot(Es, 100 * rs * np.sqrt(Es), "-", color=c, lw=1.5, label=f"sim σ_E = {g:.1f} %")
        nf = np.array([(scan, E) in NOFIT for E in ENERGIES])
        ax.errorbar(Es[~nf], 100 * rd[~nf] * np.sqrt(Es[~nf]), yerr=100 * ed[~nf] * np.sqrt(Es[~nf]), fmt="o",
                    color=ps.DATA, ms=6, capsize=3, label="data")
        if nf.any():
            ax.errorbar(Es[nf], 100 * rd[nf] * np.sqrt(Es[nf]), yerr=100 * ed[nf] * np.sqrt(Es[nf]), fmt="o",
                        mfc="none", color=ps.DATA, ms=6, capsize=3, label="data, runs ≤ 20 (not fitted)")
        ax.set_xscale("log")
        ax.set_xticks(ENERGIES)
        ax.set_xticklabels([f"{E:g}" for E in ENERGIES])
        ax.set_xlabel("beam energy [GeV]")
        ax.set_title(f"th230 {scan}", fontsize=10)
        ax.legend(fontsize=7)
    axs[0].set_ylabel("σ/μ · √E  [% √GeV]")
    fig.suptitle("Core resolution × √E: rising curve = constant term; beam spread enters it in quadrature", fontsize=10)
    fig.tight_layout()
    fig.savefig(f"{args.out}/resolution_vs_sigmaE.png", dpi=130)
    plt.close(fig)

    # ---------- figure 4: distributions at 20 and 74 GeV (E2), data vs sigma_E = 1.0 and 2.0 %
    fig, axs = plt.subplots(2, 3, figsize=(15, 6.5), gridspec_kw={"height_ratios": [3, 1.3]}, sharex="col")
    for j, E in enumerate((20, 52, 74)):
        r = R[("E2", E)]
        md, rd = core(r["data"])
        edges = np.linspace(1 - 4 * rd, 1 + 3 * rd, 50)
        ctr = 0.5 * (edges[1:] + edges[:-1])
        cd, _ = np.histogram(r["data"] / md, edges)
        nd = cd.sum()
        axs[0, j].errorbar(ctr, cd / nd, yerr=np.sqrt(cd) / nd, fmt="o", color=ps.DATA, ms=4, label="data")
        for g, c, ls in ((1.0, ps.NODIGI, "--"), (2.0, ps.DIGI, "-"), (2.6, "0.4", ":")):
            x = emulate(r["sim"], g / 100, r["dz"])
            ms_, _ = core(x)
            cs, _ = np.histogram(x / ms_, edges)
            axs[0, j].stairs(cs / cs.sum(), edges, color=c, ls=ls, lw=1.6, label=f"sim σ_E = {g:.1f} %")
            ok = cd > 0
            axs[1, j].plot(ctr[ok], (cs / cs.sum())[ok] / (cd / nd)[ok], ls=ls, color=c, lw=1.4)
        axs[0, j].axvspan(1 - 1.5 * rd, 1 + 2 * rd, color="0.94", zorder=0)
        axs[0, j].set_title(f"th230 E2  {E} GeV (grey: χ² window)", fontsize=10)
        axs[0, j].legend(fontsize=7)
        axs[1, j].axhline(1, color="k", lw=0.8)
        axs[1, j].set_ylim(0.5, 1.5)
        axs[1, j].set_xlabel("E / μ_core")
    axs[1, 0].set_ylabel("sim / data")
    axs[0, 0].set_ylabel("fraction / bin")
    fig.tight_layout()
    fig.savefig(f"{args.out}/distributions_E2.png", dpi=130)
    plt.close(fig)
    correlation(R, args.out)
    print(f"written to {args.out}")


def correlation(R, out):
    """Best sigma_E against beam energy on the fitted points: weighted constant and straight-line fits (errors scaled
    by sqrt(chi2/ndf) of the constant fit when it is > 1), Pearson and Spearman coefficients."""
    from scipy.stats import pearsonr, spearmanr
    pts = sorted(k for k in R if k not in NOFIT and np.isfinite(R[k]["best_h"][1]))
    lines = ["", f"# correlation best sigma_E vs E, {len(pts)} fitted points: " + " ".join(f"{s}:{e:g}" for s, e in pts)]
    fig, axs = plt.subplots(1, 2, figsize=(12, 4.4), sharey=True)
    for ax, key, name, col in ((axs[0], "best_h", "shape χ²", ps.DIGI), (axs[1], "best_r", "σ/μ χ²", ps.NODIGI)):
        P = [k for k in pts if np.isfinite(R[k][key][1]) and R[k][key][1] > 0]
        E = np.array([k[1] for k in P])
        y = np.array([R[k][key][0] for k in P])
        e = np.array([R[k][key][1] for k in P])
        w = 1 / e ** 2
        c0 = (w * y).sum() / w.sum()
        chi0 = float((w * (y - c0) ** 2).sum())
        scale = max(1.0, np.sqrt(chi0 / max(len(y) - 1, 1)))
        (p1, p0), cov = np.polyfit(E, y, 1, w=1 / e, cov="unscaled")
        chi1 = float((w * (y - (p0 + p1 * E)) ** 2).sum())
        sp1 = np.sqrt(cov[0, 0]) * max(1.0, np.sqrt(chi1 / max(len(y) - 2, 1)))
        r, pr = pearsonr(E, y)
        rs, prs = spearmanr(E, y)
        lines.append(f"{name}: constant {c0:.2f} +- {scale/np.sqrt(w.sum()):.2f} % (chi2/ndf {chi0:.1f}/{len(y)-1}); "
                     f"line {p0:.2f} + {100*p1:.2f}e-2 * E/GeV, slope {100*p1:.2f} +- {100*sp1:.2f} e-2 %/GeV "
                     f"({abs(p1)/sp1:.1f} sigma, chi2/ndf {chi1:.1f}/{len(y)-2}); Pearson r = {r:+.2f} (p {pr:.2f}), "
                     f"Spearman {rs:+.2f} (p {prs:.2f})")
        for (sc, En), yy, ee in zip(P, y, e):
            ax.errorbar(En, yy, yerr=ee * scale, fmt="s" if sc == "E1" else "o", color=col, ms=6, capsize=3,
                        label=f"th230 {sc}" if (sc, En) == next(k for k in P if k[0] == sc) else None)
        xx = np.linspace(5, 80, 50)
        ax.plot(xx, p0 + p1 * xx, color=col, lw=1.2, label=f"line: slope {100*p1:+.2f} ± {100*sp1:.2f} ×10⁻² %/GeV")
        ax.axhline(c0, color="k", lw=1, ls="--", label=f"constant {c0:.2f} %")
        ax.axhspan(1.0, 2.0, color="0.92", zorder=0)
        ax.set_xlabel("beam energy [GeV]")
        ax.set_title(f"{name}: Pearson r = {r:+.2f} (p = {pr:.2f}); errors ×{scale:.1f}", fontsize=10)
        ax.legend(fontsize=7)
    axs[0].set_ylabel("best beam σ_E [%]")
    axs[0].set_ylim(0, 3.5)
    fig.suptitle("Energy dependence of the fitted beam spread, runs after the collimator closure (> 20)", fontsize=10)
    fig.tight_layout()
    fig.savefig(f"{out}/sigmaE_vs_energy.png", dpi=130)
    plt.close(fig)
    open(f"{out}/beam_energy_spread.txt", "a").write("\n".join(lines) + "\n")
    print("\n".join(lines))


if __name__ == "__main__":
    main()

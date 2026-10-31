#!/usr/bin/env python
"""
nhits and E/hit, data vs simulation, against the beam energy; is the sim - data difference quadratic in E?

    python analysis/energy_trend.py <outdir> [--campaign <sim dir>] [--tag v7] [--nmax 100000] [--edge-cut 30]

Per beam point (analysis/beam_points.txt groups; th210 central scan and the th230 E1-E5 position scans), shower-like
events (nhit_chan > 0.5 x p90, as everywhere in sim_tuning), from the valtrees:
  nhits   mean nhit_chan                       (stat error sigma/sqrt(N))
  E/hit   sum(sum_energy) / sum(nhit_chan)     [MIP] (stat error by the delta method)
th230: the five positions E1-E5 are averaged per energy; their spread (RMS / sqrt(n)) is the point-to-point
(non-statistical) uncertainty, and its median relative size is also assigned to the th210 points.
Fits of D(E) = sim - data: pol1 a + bE, pol2 a + bE + cE^2, pure quadratic cE^2, and a + cE^2; chi2/ndf,
F-test pol2 vs pol1, significance of c.
--edge-cut D: keep only events whose energy barycentre (bar_x, bar_y) is more than D mm from every edge of the
active area (pad centres +-86.325 mm plus half a pad: |bar_x|, |bar_y| < 89.09 - D), after the shower-like cut.
Writes energy_trend.txt, energy_trend.npz, nhits_vs_E.png, dnhits_vs_E.png, deph_vs_E.png.
"""
import argparse
import os
import sys

import numpy as np
import uproot
from scipy import stats
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from analysis import sim_tuning as st  # noqa: E402

DATA_C, SIM_C, LIN_C, QUAD_C = "#0b0b0b", "#2a78d6", "#8c8c8c", "#eb6834"


HALF = 86.325 + 5.53 / 2          # edge of the active area [mm]


def point(paths, nmax, edge_cut=None):
    nh, se, bx, by, left = [], [], [], [], nmax
    for p in paths:
        if left <= 0:
            break
        a = uproot.open(p)["ecal"].arrays(["nhit_chan", "sum_energy", "bar_x", "bar_y"], entry_stop=left,
                                          library="np")
        left -= len(a["nhit_chan"])
        nh.append(a["nhit_chan"].astype(float)); se.append(a["sum_energy"].astype(float))
        bx.append(a["bar_x"].astype(float)); by.append(a["bar_y"].astype(float))
    nh, se, bx, by = (np.concatenate(v) for v in (nh, se, bx, by))
    k = nh > 0.5 * np.percentile(nh, 90)
    if edge_cut is not None:
        k &= (np.abs(bx) < HALF - edge_cut) & (np.abs(by) < HALF - edge_cut)
    nh, se = nh[k], se[k]
    n = len(nh)
    r = se.sum() / nh.sum()
    # delta method for the ratio of sums
    cov = np.cov(se, nh)
    r_err = abs(r) * np.sqrt(cov[0, 0] / se.mean() ** 2 + cov[1, 1] / nh.mean() ** 2
                             - 2 * cov[0, 1] / (se.mean() * nh.mean())) / np.sqrt(n)
    return {"n": n, "nh": nh.mean(), "nh_err": nh.std() / np.sqrt(n), "eph": r, "eph_err": r_err}


def fits(E, y, err):
    """Weighted fits of y(E); returns dict model -> (params, chi2, ndf)."""
    out = {}
    w = 1 / err ** 2
    for name, cols in (("pol1", [np.ones_like(E), E]), ("pol2", [np.ones_like(E), E, E ** 2]),
                       ("cE2", [E ** 2]), ("a+cE2", [np.ones_like(E), E ** 2])):
        A = np.vstack(cols).T
        Aw = A * np.sqrt(w)[:, None]
        p, *_ = np.linalg.lstsq(Aw, y * np.sqrt(w), rcond=None)
        cov = np.linalg.inv(Aw.T @ Aw)
        chi2 = float((((A @ p) - y) ** 2 * w).sum())
        out[name] = (p, np.sqrt(np.diag(cov)), chi2, len(E) - len(p))
    return out


def model(name, p, E):
    return {"pol1": lambda: p[0] + p[1] * E, "pol2": lambda: p[0] + p[1] * E + p[2] * E ** 2,
            "cE2": lambda: p[0] * E ** 2, "a+cE2": lambda: p[0] + p[1] * E ** 2}[name]()


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("outdir")
    ap.add_argument("--campaign", default=st.V7)
    ap.add_argument("--tag", default="final_v7")
    ap.add_argument("--nmax", type=int, default=100000)
    ap.add_argument("--edge-cut", type=float, default=None)
    args = ap.parse_args()
    os.makedirs(args.outdir, exist_ok=True)
    G = st.load_groups()
    sets = {"th210 (central scan)": [g for g in G if g["set"] == "th210" and
                                    g["name"].split("_")[1] in ("e7.5", "e10", "e20", "e34", "e52", "e74")
                                    or g["name"] == "th210_e99_22_-48"],
            "th230 (positions E1-E5)": [g for g in G if g["set"] == "th230" and "_E" in g["name"]
                                       and g["name"].count("_") == 2]}
    res, lines = {}, []
    for sname, gs in sets.items():
        rows = []
        for g in gs:
            sp = st.sim_paths(g["name"], args.campaign)
            dp = st.data_paths(g)
            if not sp or not dp:
                continue
            d, s = point(dp, args.nmax, args.edge_cut), point(sp, args.nmax, args.edge_cut)
            if min(d["n"], s["n"]) < 500:
                print(f"{g['name']}: skipped, {d['n']} data / {s['n']} sim events after the cuts")
                continue
            rows.append((float(g["E"]), g["name"], d, s))
            print(f"{g['name']:<22} E {g['E']:>5}  nhits data {d['nh']:7.1f} sim {s['nh']:7.1f}  "
                  f"E/hit data {d['eph']:.3f} sim {s['eph']:.3f}  (N {d['n']}, {s['n']})", flush=True)
        Es = sorted({r[0] for r in rows})
        P = {k: [] for k in ("E", "nh_d", "nh_s", "dnh", "dnh_stat", "dnh_pp", "eph_d", "eph_s", "deph",
                             "deph_stat", "deph_pp", "npos")}
        for E in Es:
            rr = [r for r in rows if r[0] == E]
            dnh = np.array([r[3]["nh"] - r[2]["nh"] for r in rr])
            deph = np.array([r[3]["eph"] - r[2]["eph"] for r in rr])
            P["E"].append(E); P["npos"].append(len(rr))
            P["nh_d"].append(np.mean([r[2]["nh"] for r in rr])); P["nh_s"].append(np.mean([r[3]["nh"] for r in rr]))
            P["eph_d"].append(np.mean([r[2]["eph"] for r in rr])); P["eph_s"].append(np.mean([r[3]["eph"] for r in rr]))
            P["dnh"].append(dnh.mean()); P["deph"].append(deph.mean())
            P["dnh_stat"].append(np.sqrt(np.mean([r[2]["nh_err"] ** 2 + r[3]["nh_err"] ** 2 for r in rr]) / len(rr)))
            P["deph_stat"].append(np.sqrt(np.mean([r[2]["eph_err"] ** 2 + r[3]["eph_err"] ** 2 for r in rr]) / len(rr)))
            P["dnh_pp"].append(dnh.std(ddof=1) / np.sqrt(len(rr)) if len(rr) > 1 else np.nan)
            P["deph_pp"].append(deph.std(ddof=1) / np.sqrt(len(rr)) if len(rr) > 1 else np.nan)
        res[sname] = {k: np.array(v) for k, v in P.items()}
    # point-to-point uncertainty: relative size from th230 (spread over positions / nhits), applied where missing
    t = res["th230 (positions E1-E5)"]
    rel_nh = np.nanmedian(t["dnh_pp"] * np.sqrt(t["npos"]) / t["nh_d"])        # per single position
    rel_eph = np.nanmedian(t["deph_pp"] * np.sqrt(t["npos"]) / t["eph_d"])
    lines.append(f"point-to-point (position-to-position) spread of a single point: nhits {rel_nh * 100:.2f} %, "
                 f"E/hit {rel_eph * 100:.2f} % (median over th230 energies)")
    for sname, r in res.items():
        for q, rel in (("dnh", rel_nh), ("deph", rel_eph)):
            base = r["nh_d"] if q == "dnh" else r["eph_d"]
            pp = np.where(np.isfinite(r[f"{q}_pp"]), r[f"{q}_pp"], rel * base / np.sqrt(r["npos"]))
            r[f"{q}_err"] = np.sqrt(r[f"{q}_stat"] ** 2 + pp ** 2)
            F = fits(r["E"], r[q], r[f"{q}_err"])
            r[f"{q}_fits"] = F
            c1, n1 = F["pol1"][2], F["pol1"][3]
            c2, n2 = F["pol2"][2], F["pol2"][3]
            fstat = (c1 - c2) / (c2 / n2) if n2 > 0 and c2 > 0 else np.nan
            pval = 1 - stats.f.cdf(fstat, 1, n2) if np.isfinite(fstat) else np.nan
            lab = "nhits" if q == "dnh" else "E/hit [MIP]"
            lines.append(f"\n{sname}  D {lab} = sim - data   E: {' '.join(f'{e:g}' for e in r['E'])}")
            lines.append("  D: " + " ".join(f"{v:+.3f}" if q == "deph" else f"{v:+.1f}" for v in r[q])
                         + "   err: " + " ".join(f"{v:.3f}" if q == "deph" else f"{v:.1f}" for v in r[f"{q}_err"]))
            for m, (p, pe, c, n) in F.items():
                lines.append(f"  {m:<6} chi2/ndf {c:8.1f}/{n}  p-value {1 - stats.chi2.cdf(c, n) if n > 0 else np.nan:.3g}"
                             f"  params " + " ".join(f"{a:+.4g}+-{b:.2g}" for a, b in zip(p, pe)))
            lines.append(f"  pol2 vs pol1: F = {fstat:.2f}, p = {pval:.3g};  c/sigma_c = "
                         f"{F['pol2'][0][2] / F['pol2'][1][2]:+.1f}")
    txt = "\n".join(lines)
    print(txt)
    open(f"{args.outdir}/energy_trend.txt", "w").write(txt + "\n")
    np.savez(f"{args.outdir}/energy_trend.npz", **{f"{s}|{k}": v for s, r in res.items() for k, v in r.items()
                                                    if not k.endswith("_fits")})

    plt.rcParams.update({"axes.spines.top": False, "axes.spines.right": False, "axes.grid": True,
                         "grid.color": "#e6e6e3", "grid.linewidth": 0.8, "axes.edgecolor": "#8c8c8c",
                         "axes.labelcolor": "#0b0b0b", "font.size": 10})
    Ef = np.linspace(0, 105, 200)
    sel = "" if args.edge_cut is None else f", barycentre > {args.edge_cut:g} mm from the edges"
    # 1. nhits vs E
    fig, ax = plt.subplots(1, 2, figsize=(11, 4.2))
    for a, (sname, r) in zip(ax, res.items()):
        a.errorbar(r["E"], r["nh_d"], yerr=r["dnh_err"] / np.sqrt(2), fmt="o", ms=6, color=DATA_C, label="data", zorder=3)
        a.plot(r["E"], r["nh_s"], "-s", ms=6, lw=2, color=SIM_C, label=f"simulation ({args.tag})", zorder=2)
        a.set_title(sname + sel, fontsize=10, loc="left"); a.set_xlabel("beam energy [GeV]"); a.set_xlim(0, 105)
    ax[0].set_ylabel("hits per event (mean, shower-like events)"); ax[0].legend(frameon=False)
    fig.tight_layout(); fig.savefig(f"{args.outdir}/nhits_vs_E.png", dpi=130); plt.close(fig)
    # 2./3. differences with fits
    for q, fn, ylab in (("dnh", "dnhits_vs_E.png", "nhits(sim) - nhits(data)"),
                        ("deph", "deph_vs_E.png", "E/hit(sim) - E/hit(data)  [MIP]")):
        fig, ax = plt.subplots(1, 2, figsize=(11, 4.2))
        for a, (sname, r) in zip(ax, res.items()):
            F = r[f"{q}_fits"]
            a.axhline(0, color="#52514e", lw=1)
            for m, col, ls in (("pol1", LIN_C, "--"), ("pol2", QUAD_C, "-")):
                p, pe, c, n = F[m]
                a.plot(Ef, model(m, p, Ef), ls, lw=2, color=col,
                       label=f"{'linear' if m == 'pol1' else 'quadratic'}: chi2/ndf {c:.1f}/{n}")
            a.errorbar(r["E"], r[q], yerr=r[f"{q}_err"], fmt="o", ms=6, color=DATA_C, zorder=3,
                       label="sim - data")
            a.set_title(sname + sel, fontsize=10, loc="left"); a.set_xlabel("beam energy [GeV]"); a.set_xlim(0, 105)
            a.legend(frameon=False, fontsize=8)
        ax[0].set_ylabel(ylab)
        fig.tight_layout(); fig.savefig(f"{args.outdir}/{fn}", dpi=130); plt.close(fig)


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Energy linearity and resolution of one energy scan at a fixed table position, data against both simulations,
from the per-group summary of final_panels.py (Gaussian-core mu and sigma/mu of the calibrated energy).

The resolution is only meaningful where the showers stay away from dead regions: the energy scans at E1 and E2
(th230, runs 12-29), and the th210 scan S1 (eudaq 285-293, the only normal-incidence th210 scan; no th210 run was
taken at E1 or E2).

    python3 analysis/position_scans.py <final_panels outdir> [--scan th230:E1 --scan th230:E2 --scan th210:S1]
Writes lin_<set>_<pos>.png, res_<set>_<pos>.png and position_scans.txt in --out (default: the input directory).
Statistical errors: sigma/sqrt(N) on mu, (sigma/mu)/sqrt(2N) on sigma/mu.

Calibration: by default the energies are final_panels.py's, i.e. the MIP -> GeV line of each sample fitted on the th210
scan. --recal <set>:<pos> --recal-dir <no-cut final_panels outdir> refits that straight line, per sample, on the given
scan without the border cut (same method as final_panels.py: unweighted straight line through the Gaussian-core mu of
sum_w against E) and recalibrates every point with it; the output files get the suffix _cal<set><pos>. sum_w is
recovered from the th210-calibrated values through the th210 line (<recal-dir>/mip_gev_line.txt); the resolution is
recomputed with the new line too (sigma/mu of the energy depends on the line's offset).
"""
import argparse
import os
import re
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from scipy.optimize import curve_fit  # noqa: E402

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from analysis import plot_style as ps  # noqa: E402

KINDS = ("data", "nodigi", "digi")
LABEL = {"data": "Data", "nodigi": "Sim (no digitisation)", "digi": "Sim (digitised)"}
COL = {"data": ps.DATA, "nodigi": ps.NODIGI, "digi": ps.DIGI}
LS = {"data": "-", "nodigi": ":", "digi": "--"}
# group-name pattern of each scan in beam_points.txt
MS, CAP = 2.5, 2.5   # small markers so the (small) statistical error bars and their caps stay visible


def data_errors(ax, x, y, yerr):
    """Data error bars with caps, drawn on their own: data_points uses mew=0, which also hides the caps."""
    ax.errorbar(x, y, yerr=yerr, fmt="none", ecolor=ps.DATA, elinewidth=1.2, capsize=CAP, capthick=1.2, zorder=4)
SCANS = {"S1": r"^th210_e(7\.5|10|20|34|52|74|99)_\d+_-\d+$"}


def read_summary(path):
    out, g = {}, None
    for line in open(path):
        m = re.match(r"== (\S+)\s+\(.*N data/nodigi/digi = (\d+)/(\d+)/(\d+)", line)
        if m:
            g = m.group(1)
            out[g] = {"N": tuple(int(m.group(i)) for i in (2, 3, 4))}
            continue
        f = line.split()
        if g and f and f[0] == "gev":
            out[g]["mu"] = tuple(float(v) for v in f[1:4])
            out[g]["res"] = tuple(float(v.rstrip("%")) / 100 for v in f[4:7])
    return out


def read_points(path):
    pts = {}
    for line in open(path):
        f = line.split("#")[0].split()
        if len(f) >= 6:
            pts[f[0]] = dict(set=f[1], pos=f[2], E=float(f[3]), x=float(f[4]), y=float(f[5]), runs=f[8] if len(f) > 8 else "")
    return pts


def stoch(E, a, c):
    return np.hypot(a / np.sqrt(E), c)


def scan_groups(pts, summ, st, pos):
    if pos in SCANS:
        names = [g for g in pts if pts[g]["set"] == st and re.match(SCANS[pos], g) and pts[g]["y"] < -45]
    else:
        names = [g for g in pts if pts[g]["set"] == st and pts[g]["pos"] == pos]
    names = [g for g in names if g in summ and "mu" in summ[g]]
    return sorted(names, key=lambda g: pts[g]["E"])


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("indir")
    ap.add_argument("--scan", action="append", default=None, help="<set>:<position>, repeatable")
    ap.add_argument("--points", default=os.path.join(os.path.dirname(os.path.abspath(__file__)), "beam_points.txt"))
    ap.add_argument("--title", default="", help="selection note for the figure titles")
    ap.add_argument("--recal", default=None, help="<set>:<pos> scan to refit the MIP -> GeV line on (e.g. th230:E2)")
    ap.add_argument("--recal-dir", default=None, help="final_panels outdir without the border cut (for --recal)")
    ap.add_argument("--out", default=None)
    args = ap.parse_args(argv)
    scans = args.scan or ["th230:E1", "th230:E2"]
    out = args.out or args.indir
    pts, summ = read_points(args.points), read_summary(f"{args.indir}/final_panels_summary.txt")
    suffix = ""
    if args.recal:
        rst, rpos = args.recal.split(":")
        base = read_summary(f"{args.recal_dir}/final_panels_summary.txt")
        line210 = {}
        for l in open(f"{args.recal_dir}/mip_gev_line.txt"):
            f = l.split()
            if f and not l.startswith("#"):
                line210[f[0]] = (float(f[1]), float(f[2]))
        rnames = scan_groups(pts, base, rst, rpos)
        newline = {}
        for i, k in enumerate(KINDS):
            a0, b0 = line210[k]
            Er = np.array([pts[g]["E"] for g in rnames])
            mw = np.array([a0 * base[g]["mu"][i] + b0 for g in rnames])
            newline[k] = np.polyfit(Er, mw, 1)
        for g, s_ in summ.items():
            if "mu" not in s_:
                continue
            mu, res = [], []
            for i, k in enumerate(KINDS):
                a0, b0 = line210[k]
                a1, b1 = newline[k]
                mw, sw = a0 * s_["mu"][i] + b0, a0 * s_["mu"][i] * s_["res"][i]
                m1 = (mw - b1) / a1
                mu.append(m1)
                res.append(sw / a1 / m1)
            s_["mu"], s_["res"] = tuple(mu), tuple(res)
        suffix = f"_cal{rst}{rpos}"
        note = (f"MIP -> GeV line refitted per sample on {rst} {rpos} without the border cut ({', '.join(rnames)}): "
                + "; ".join(f"{k} a = {newline[k][0]:.2f} MIP/GeV, b = {newline[k][1]:.1f} MIP" for k in KINDS)
                + "; th210 line: " + "; ".join(f"{k} a = {line210[k][0]:.2f}, b = {line210[k][1]:.1f}" for k in KINDS))
        print(note)
        args.title = (args.title + "\n" if args.title else "") + f"MIP→GeV line of each sample fitted on {rst} {rpos} (no border cut)"
    plt.rcParams.update(ps.RC)
    lines = [f"# energy scans at a fixed position; {args.title}",
             "# scan group E N_data | mu/E data nodigi digi | sigma/mu data nodigi digi"]
    if args.recal:
        lines.insert(1, "# " + note)
    Ef = np.linspace(5, 105, 200)
    for sc in scans:
        st, pos = sc.split(":")
        names = scan_groups(pts, summ, st, pos)
        if not names:
            print(f"{sc}: no groups")
            continue
        E = np.array([pts[g]["E"] for g in names])
        runs = ", ".join(pts[g]["runs"] for g in names)
        val = {}
        for i, k in enumerate(KINDS):
            N = np.array([summ[g]["N"][i] for g in names], float)
            mu = np.array([summ[g]["mu"][i] for g in names])
            res = np.array([summ[g]["res"][i] for g in names])
            val[k] = (mu, mu * res / np.sqrt(N), res, res / np.sqrt(2 * N))
        for g, e in zip(names, E):
            s = summ[g]
            lines.append(f"{sc} {g} {e:g} {s['N'][0]} | " + " ".join(f"{m / e:.3f}" for m in s["mu"]) + " | "
                         + " ".join(f"{100 * r:.2f}%" for r in s["res"]))
        tag = f"{st}_{pos}"
        where = f"{st}, position {pos}" + (" (th210 energy scan, eudaq 285-293)" if pos == "S1" else "")

        # linearity: mu against E (top), (mu - E)/E (bottom)
        fig, (a1, a2) = plt.subplots(2, 1, figsize=(6.4, 5.6), sharex=True, gridspec_kw=dict(height_ratios=[2, 1.2]))
        for k in KINDS:
            mu, dmu, _, _ = val[k]
            if k == "data":
                ps.data_points(a1, E, mu, yerr=dmu, label=LABEL[k] + " (stat. errors)", ms=MS)
                data_errors(a1, E, mu, dmu)
                ps.data_points(a2, E, mu / E - 1, yerr=dmu / E, label=None, ms=MS)
                data_errors(a2, E, mu / E - 1, dmu / E)
            else:
                ps.sim_line(a1, E, mu, dict(color=COL[k], ls=LS[k], marker=None), label=LABEL[k])
                ps.sim_line(a2, E, mu / E - 1, dict(color=COL[k], ls=LS[k], marker=None))
                a2.errorbar(E, mu / E - 1, yerr=dmu / E, color=COL[k], ls="none", marker="o" if k == "digi" else "s",
                            ms=MS, mfc=COL[k] if k == "digi" else "white", capsize=CAP, elinewidth=1)
        a1.plot([0, 105], [0, 105], color="0.5", lw=0.8, label="$\\mu$ = E")
        a1.set_ylabel("calibrated energy $\\mu$ [GeV]")
        a1.legend(frameon=False, fontsize=8)
        a2.text(0.99, 0.04, "error bars: statistical only", transform=a2.transAxes, ha="right", fontsize=7.5, color="0.35")
        a2.axhline(0, color="0.5", lw=0.8)
        a2.set_ylabel("($\\mu$ − E) / E")
        a2.set_xlabel("beam energy E [GeV]")
        lo = min(-0.005, np.nanmin([val[k][0] / E - 1 for k in KINDS]) - 0.008)
        hi = max(0.005, np.nanmax([val[k][0] / E - 1 for k in KINDS]) + 0.008)
        a2.set_ylim(lo, hi)
        for a in (a1, a2):
            a.grid(alpha=0.25)
        a1.set_title(f"Linearity, {where}\n{args.title}", fontsize=9)
        fig.tight_layout()
        fig.savefig(f"{out}/lin_{tag}{suffix}.png", dpi=130)
        plt.close(fig)

        # resolution: sigma/mu against E with a/sqrt(E) (+) c (top), sim/data (bottom)
        fig, (a1, a2) = plt.subplots(2, 1, figsize=(6.4, 5.6), sharex=True, gridspec_kw=dict(height_ratios=[2, 1.2]))
        for k in KINDS:
            _, _, res, dres = val[k]
            try:
                (a, c), _ = curve_fit(stoch, E, res, sigma=dres, p0=(0.2, 0.01))
                fl = f": {100 * abs(a):.1f} %/$\\sqrt{{E}}$ $\\oplus$ {100 * abs(c):.2f} %"
                lines.append(f"fit {sc} {k:7s} a = {100 * abs(a):5.1f} %/sqrt(GeV)  c = {100 * abs(c):4.2f} %")
                a1.plot(Ef, 100 * stoch(Ef, a, c), color=COL[k], lw=1, ls=LS[k], alpha=0.8)
            except RuntimeError:
                fl = ""
            if k == "data":
                ps.data_points(a1, E, 100 * res, yerr=100 * dres, label=LABEL[k] + fl, ms=MS)
                data_errors(a1, E, 100 * res, 100 * dres)
            else:
                a1.errorbar(E, 100 * res, yerr=100 * dres, color=COL[k], ls="none", capsize=CAP, elinewidth=1,
                            marker="o" if k == "digi" else "s", ms=MS, mfc=COL[k] if k == "digi" else "white",
                            label=LABEL[k] + fl)
                r = res / val["data"][2]
                dr = r * np.hypot(dres / res, val["data"][3] / val["data"][2])
                ps.sim_line(a2, E, r, dict(color=COL[k], ls=LS[k], marker=None))
                a2.errorbar(E, r, yerr=dr, color=COL[k], ls="none", marker="o" if k == "digi" else "s", ms=MS,
                            mfc=COL[k] if k == "digi" else "white", capsize=CAP, elinewidth=1)
        a1.set_ylabel("energy resolution $\\sigma/\\mu$ [%]")
        a1.set_ylim(0, None)
        a1.legend(frameon=False, fontsize=8, title="error bars: statistical only", title_fontsize=7.5)
        a2.axhline(1, color="0.5", lw=0.8)
        a2.set_ylabel("sim / data")
        a2.set_xlabel("beam energy E [GeV]")
        for a in (a1, a2):
            a.grid(alpha=0.25)
        a1.set_title(f"Energy resolution, {where}\n{args.title}", fontsize=9)
        fig.tight_layout()
        fig.savefig(f"{out}/res_{tag}{suffix}.png", dpi=130)
        plt.close(fig)
        print(f"{sc}: {len(names)} groups, runs {runs}")
    open(f"{out}/position_scans{suffix}.txt", "w").write("\n".join(lines) + "\n")
    print("\n".join(lines))


if __name__ == "__main__":
    sys.exit(main())

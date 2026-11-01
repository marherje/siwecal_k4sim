#!/usr/bin/env python3
"""Energy resolution and linearity of the laterally contained groups, from the per-group summaries of
final_panels.py (with and without --border-cut).

A group enters when its beam centre (beam_points.txt) is inside the cut, min(L - |x|, L - |y|) >= cut - margin, and the
cut leaves at least --min-data data events. For the other groups the cut keeps only the tail of the beam spot (and,
at high energy where the spot is narrow, mostly particles that are not beam electrons), so their numbers do not
describe the shower. Fits sigma/E = a/sqrt(E) (+) c per sample on the contained groups.

    python3 analysis/contained_resolution.py <b50 dir> <nocut dir> --cut 50 [--out <dir>]
Writes contained_resolution.png and contained_resolution.txt.
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

L_ACTIVE = 89.09   # half-size of the active area [mm], as barxy_min_border_distance
KINDS = ("data", "nodigi", "digi")
LABEL = {"data": "Data", "nodigi": "Sim (no digitisation)", "digi": "Sim (digitised)"}
MARK = {"th210": "o", "th220": "s", "th230": "^"}


def read_points(path):
    pts = {}
    for line in open(path):
        f = line.split("#")[0].split()
        if len(f) >= 6:
            pts[f[0]] = dict(set=f[1], E=float(f[3]), x=float(f[4]), y=float(f[5]))
    return pts


def read_summary(path):
    """{group: {'N': (nd, nn, ng), 'gev': (mu_d, mu_n, mu_g, res_d, res_n, res_g)}} from final_panels_summary.txt."""
    out, g = {}, None
    for line in open(path):
        m = re.match(r"== (\S+)\s+\(.*N data/nodigi/digi = (\d+)/(\d+)/(\d+)", line)
        if m:
            g = m.group(1)
            out[g] = {"N": tuple(int(m.group(i)) for i in (2, 3, 4))}
            continue
        f = line.split()
        if g and f and f[0] == "gev":
            out[g]["gev"] = tuple(float(v) for v in f[1:4]) + tuple(float(v.rstrip("%")) / 100 for v in f[4:7])
    return out


def stoch(E, a, c):
    return np.hypot(a / np.sqrt(E), c)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("cut_dir")
    ap.add_argument("nocut_dir")
    ap.add_argument("--cut", type=float, default=50.0)
    ap.add_argument("--min-data", type=int, default=1000)
    ap.add_argument("--margin", type=float, default=0.0,
                    help="accept beam centres up to this far outside the cut [mm] (cut - margin <= centre distance)")
    ap.add_argument("--points", default=os.path.join(os.path.dirname(os.path.abspath(__file__)), "beam_points.txt"))
    ap.add_argument("--uniformity", action="store_true", help="also the th220 uniformity scan (U2-U50, all 74 GeV)")
    ap.add_argument("--outlier", type=float, default=1.3,
                    help="groups whose data sigma/mu exceeds this times the digitised one are drawn open and left out "
                         "of the fits (data-side problems: dead regions, hadrons at 99 GeV)")
    ap.add_argument("--out", default=None)
    args = ap.parse_args(argv)
    out = args.out or args.cut_dir
    pts = read_points(args.points)
    cut, nocut = read_summary(f"{args.cut_dir}/final_panels_summary.txt"), \
        read_summary(f"{args.nocut_dir}/final_panels_summary.txt")

    rows, lines = [], [f"# contained groups: beam centre >= {args.cut - args.margin:g} mm from the active edge (L = {L_ACTIVE} mm) "
                       f"and >= {args.min_data} data events after barxy_min_border_distance > {args.cut:g} mm",
                       "# group set E centre_dist N_data | res data nodigi digi (cut) | mu/E data nodigi digi (cut) | "
                       "res data digi (no cut)"]
    for g, s in cut.items():
        p = pts.get(g)
        if p is None or "gev" not in s or (re.search(r"_U\d+_e", g) and not args.uniformity):
            continue
        dist = min(L_ACTIVE - abs(p["x"]), L_ACTIVE - abs(p["y"]))
        if dist < args.cut - args.margin or s["N"][0] < args.min_data:
            continue
        mu, res = s["gev"][:3], s["gev"][3:]
        nc = nocut[g]["gev"]
        out_ = res[0] > args.outlier * res[2]
        rows.append((g, p["set"], p["E"], mu, res, nc, out_))
        lines.append(f"{g:24s} {p['set']} {p['E']:5g} {dist:5.1f} {s['N'][0]:7d} | "
                     + " ".join(f"{100 * r:5.2f}%" for r in res) + " | "
                     + " ".join(f"{m / p['E']:.3f}" for m in mu) + " | "
                     + f"{100 * nc[3]:5.2f}% {100 * nc[5]:5.2f}%" + ("  OUTLIER (not fitted)" if out_ else ""))

    plt.rcParams.update(ps.RC)
    sets = sorted({r[1] for r in rows})
    fig, ax = plt.subplots(1, len(sets) + 1, figsize=(5 * (len(sets) + 1), 4.3))
    col = {"data": ps.DATA, "nodigi": ps.NODIGI, "digi": ps.DIGI}
    lstyle = {"data": "-", "nodigi": ":", "digi": "--"}
    lines.append("# fit sigma/E = a/sqrt(E) (+) c on the contained groups, per sample and threshold set (outliers excluded)")
    Ef = np.linspace(6, 105, 200)
    lin = ax[-1]
    for j, st in enumerate(sets):
        a_res = ax[j]
        sel = [r for r in rows if r[1] == st]
        E = np.array([r[2] for r in sel])
        good = np.array([not r[6] for r in sel])
        ymax = 0
        for i, k in enumerate(KINDS):
            res = np.array([r[4][i] for r in sel])
            mu = np.array([r[3][i] for r in sel]) / E
            kw = dict(color=col[k], marker=MARK[st], ls="none", ms=6, mfc=col[k] if k != "nodigi" else "white", mew=1.3)
            a_res.plot(E[good], 100 * res[good], **kw)
            lin.plot(E[good], mu[good], **kw)
            ymax = max(ymax, 100 * res[good].max())
            if len(set(E[good])) >= 3:
                try:
                    (a, c), _ = curve_fit(stoch, E[good], res[good], p0=(0.2, 0.01))
                    a_res.plot(Ef, 100 * stoch(Ef, a, c), color=col[k], lw=1.1, ls=lstyle[k],
                               label=f"{LABEL[k]}: {100 * abs(a):.1f} %/$\\sqrt{{E}}$ $\\oplus$ {100 * abs(c):.2f} %")
                    lines.append(f"fit {st} {k:7s} a = {100 * abs(a):5.1f} %/sqrt(GeV)  c = {100 * abs(c):4.2f} %"
                                 f"  ({good.sum()} groups)")
                except RuntimeError:
                    pass
        a_res.set_ylim(0, ymax * 1.15)
        outl = [r for r in sel if r[6]]
        if outl:
            a_res.text(0.03, 0.04, "not fitted, data σ/μ > " + f"{args.outlier:g}" + " × digitised:\n" + "\n".join(
                f"{r[0].split('_e', 1)[1]}: {100 * r[4][0]:.1f} %" for r in outl), transform=a_res.transAxes,
                fontsize=7.5, color=ps.DATA, va="bottom")
        a_res.set_title(st)
        a_res.set_ylabel(r"energy resolution $\sigma/\mu$ [%]")
        a_res.legend(fontsize=7.5, loc="upper right", bbox_to_anchor=(1, 0.93), frameon=False)
    lin.set_ylabel(r"$\mu$ / beam energy")
    lin.axhline(1, color="0.5", lw=0.8)
    lin.set_title("linearity (outliers not drawn)")
    for a in ax:
        a.set_xlabel("beam energy [GeV]")
        a.grid(alpha=0.25)
    from matplotlib.lines import Line2D
    h = [Line2D([], [], color=col[k], marker="o", ls="none", ms=6, mfc=col[k] if k != "nodigi" else "white",
                mew=1.3, label=LABEL[k]) for k in KINDS]
    h += [Line2D([], [], color="0.35", marker=MARK[s_], ls="none", ms=6, label=s_) for s_ in sets]
    fig.legend(handles=h, loc="upper center", ncol=len(h), frameon=False, fontsize=8.5)
    fig.tight_layout(rect=[0, 0, 1, 0.92])
    fig.savefig(f"{out}/contained_resolution.png", dpi=130)
    open(f"{out}/contained_resolution.txt", "w").write("\n".join(lines) + "\n")
    print("\n".join(lines))


if __name__ == "__main__":
    sys.exit(main())

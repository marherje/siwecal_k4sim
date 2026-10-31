#!/usr/bin/env python
"""
Figures of the "beam-line material and chip occupancy" note (2026-10-26): what went into the new default
SND_compact_beamline.xml and into AdcDigitizer.OccupancySlope, and why.

    python analysis/beamline_occupancy_figs.py <outdir> [--test-dir <dir with t2/ t3/>] [--only s1,s2]

s1  material: longitudinal profile and lateral width, near vs far (sim_tuning caches of the position scan)
    profile_long.png   normalised longitudinal profile, data / v7 / 0.2 X0 Al at -2.2 m / at -30 m / H2-like
    lateral.png        radial energy density of layer 1 and R90 per layer (sim - data)
    lever.png          L1-2 widening against the distance of the material, with the multiple-scattering estimate
    delay.png          arrival time against the beam origin and the shaping window; trigger time far vs near
s2  occupancy: the measured baseline sag, the k scan, and the in-digitiser vs post-hoc check
    pedestal.png       <B> against the chip charge Q from beam data (cm_trigger_*.txt), and the turn-on shift
    kscan_<set>.png    per threshold set: chip-multiplicity chi2 and hits sim/data against k (occupancy_fit.txt)
    impl_check.png     hit survival against Q: AdcDigitizer (OccupancySlope) vs sim_tuning.apply_baseline

Style (analysis/plot_style.py): data blue points with error bars, the digitised simulations as histograms in a gradient
of reds with their own hatch; ratios and trends draw the simulations as lines and the data as points.
--test-dir needs t2/d_{0,on,nooff}/ecal.root (74 GeV e-, H2-like line) and t3/d_0/ecal.root (final_v7 geometry):
analysis/note_test_samples.sh makes them.
"""
import argparse
import os
import re
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from analysis import sim_tuning as st  # noqa: E402
from analysis.tuning_lateral import lateral, r90, RB  # noqa: E402
from analysis import plot_style as ps  # noqa: E402
from analysis.plot_style import INK, INK2, GRID  # noqa: E402

NL = st.NL
VARIANTS = ("v7", "Al -2.2 m", "Al -30 m", "H2-like")
STYLE = dict(zip(VARIANTS, ps.reds(len(VARIANTS))))     # digitised models, light (v7) to dark (H2-like)
LAB = {"data": "data", "v7": "sim_v7 (0.04 X0, last 2.6 m)", "Al -2.2 m": "+0.20 X0 Al at -2.2 m (near)",
       "Al -30 m": "+0.20 X0 Al at -30 m (far)", "H2-like": "H2-like line (0.19 X0, mean 39 m)"}
CACHE = {"v7": "pos_v7", "Al -2.2 m": "pos_Al__2_2_m", "Al -30 m": "pos_Al__30_m", "H2-like": "pos_H2_like"}
DIST = {"Al -2.2 m": 2.2, "Al -10 m": 10.0, "split": None, "Al -30 m": 30.0, "H2-like": 39.0}
CM_RUNS = {"e291": "20 GeV", "e287": "52 GeV", "e286": "74 GeV"}

plt.rcParams.update(ps.RC)


def sample(name, group):
    p = f"{st.OUT}/cache/{'data' if name == 'data' else CACHE[name]}_{group}.npz"
    return dict(np.load(p))


def long_profile(s):
    nev = int(s["nev"])
    epl = np.zeros((nev, NL))
    np.add.at(epl, (s["ev"], s["slab"].astype(int)), s["e"].astype(float))
    tot = epl.sum(1)
    f = epl[tot > 0] / tot[tot > 0, None]
    return f.mean(0), f.std(0) / np.sqrt(len(f))


# ------------------------------------------------------------------------------------------------------------- s1
def fig_profile_long(out):
    groups = [("th230_E4_e10", "10 GeV"), ("th230_E4_e20", "20 GeV"), ("th230_E4_e74", "74 GeV")]
    fig, ax = plt.subplots(2, 3, figsize=(12, 6), sharex=True, gridspec_kw={"height_ratios": [2, 1]})
    L = np.arange(NL)
    edges = np.arange(NL + 1) - 0.5
    for j, (g, lab) in enumerate(groups):
        fd, ed = long_profile(sample("data", g))
        for n in VARIANTS:
            f, _ = long_profile(sample(n, g))
            # several overlapping models: only the chosen default keeps its hatch, the others are outlines
            ps.sim_hist(ax[0, j], f, edges, STYLE[n], label=LAB[n], fill=(n == "H2-like"))
            ps.sim_hist(ax[1, j], f / fd, edges, STYLE[n], fill=False)
        ps.data_points(ax[0, j], L, fd, ed, edges, label=LAB["data"])
        ax[1, j].axhline(1, color=INK2, lw=1)
        ax[1, j].set_ylim(0.8, 1.25)
        ax[0, j].set_title(f"th230 E4, e$^-$ {lab}", fontsize=10)
        ax[1, j].set_xlabel("layer")
    ax[0, 0].set_ylabel("fraction of the event energy")
    ax[1, 0].set_ylabel("sim / data")
    fig.suptitle("Longitudinal profile: the AMOUNT matters (v7 starts too late), the distance does not "
                 "(near and far overlap)", fontsize=10)
    h, l = ax[0, 0].get_legend_handles_labels()   # one row under the title: inside a panel it hides the peak
    fig.legend(h, l, loc="upper center", ncol=len(l), fontsize=8, bbox_to_anchor=(0.5, 0.955))
    fig.tight_layout(rect=[0, 0, 1, 0.93])
    fig.savefig(f"{out}/profile_long.png", dpi=120)
    plt.close(fig)


def fig_lateral(out):
    groups = [("th230_E4_e10", "10 GeV"), ("th230_E4_e20", "20 GeV"), ("th230_E4_e74", "74 GeV")]
    fig, ax = plt.subplots(2, 3, figsize=(12, 6.2))
    rc = 0.5 * (RB[1:] + RB[:-1])
    area = np.pi * (RB[1:] ** 2 - RB[:-1] ** 2)
    ledges = np.arange(10) - 0.5
    for j, (g, lab) in enumerate(groups):
        prof = {n: lateral(sample(n, g))["prof"] for n in ("data",) + VARIANTS}
        for n in VARIANTS:
            p = prof[n]
            ps.sim_hist(ax[0, j], p[1] / p[1].sum() / area, RB, STYLE[n], label=LAB[n], fill=False)  # log axis
            R = [r90(p[l]) - r90(prof["data"][l]) for l in range(9)]
            ps.sim_hist(ax[1, j], np.array(R), ledges, STYLE[n], fill=False)
        p = prof["data"]
        ps.data_points(ax[0, j], rc, p[1] / p[1].sum() / area, None, RB, label=LAB["data"])
        ax[1, j].axhline(0, color=INK2, lw=1)
        ps.pad_band(ax[1, j], label=(j == 0))
        ax[0, j].set_yscale("log")
        ax[0, j].set_xlim(0, 80)
        ax[0, j].set_xlabel("r from the shower axis [mm]")
        ax[0, j].set_title(f"layer 1, e$^-$ {lab}", fontsize=10)
        ax[1, j].axhline(0, color=INK2, lw=0.8)
        ax[1, j].set_xlabel("layer")
    ax[0, 0].set_ylabel("energy density [fraction / mm$^2$]")
    ax[1, 0].set_ylabel("R90 sim $-$ data [mm]")
    ax[1, 0].legend(fontsize=7.5, loc="lower right")
    fig.suptitle("Lateral width of the first layers: the DISTANCE matters (only far material gives the data's halo)",
                 fontsize=10)
    h, l = ax[0, 0].get_legend_handles_labels()
    fig.legend(h, l, loc="upper center", ncol=len(l), fontsize=8, bbox_to_anchor=(0.5, 0.955))
    fig.tight_layout(rect=[0, 0, 1, 0.93])
    fig.savefig(f"{out}/lateral.png", dpi=120)
    plt.close(fig)


def fig_lever(out):
    """Measured L1-2 R90 widening (variant - v7) vs distance, with theta0 * L for 0.2 X0."""
    rows = {}
    for line in open(f"{st.OUT}/upstream/position/position.txt"):
        m = re.match(r"(\S+) (.+?)\s+([\d.]+)\s+([\d.]+)\s+(\d+) \| (.+?) \|", line)
        if m:
            rows[(m.group(1), m.group(2).strip())] = np.array(m.group(6).split(), float)
    fig, ax = plt.subplots(1, 2, figsize=(11, 4.2))
    zz = np.linspace(0, 45, 200)
    energies = (("th230_E4_e10", 10), ("th230_E4_e20", 20), ("th230_E4_e74", 74))
    for (g, E), sty, mk in zip(energies, ps.reds(3), ("o", "s", "^")):
        v7 = rows[(g, "v7")][1:3].mean()
        xs, ys = [], []
        for n, d in DIST.items():
            if d is not None and (g, n) in rows:
                xs.append(d); ys.append(rows[(g, n)][1:3].mean() - v7)
        o = np.argsort(xs)
        ps.sim_line(ax[0], np.array(xs)[o], np.array(ys)[o], dict(sty, marker=mk), label=f"sim, {E} GeV")
        need = rows[(g, "data")][1:3].mean() - v7
        ax[0].axhline(need, color=ps.DATA, lw=1.1, ls=sty["ls"])
        ax[0].text(44, need, f"data {E} GeV", color=ps.DATA, fontsize=7.5, ha="right", va="bottom")
        th0 = 13.6e-3 / E * np.sqrt(0.2) * (1 + 0.038 * np.log(0.2))      # rad
        ax[1].plot(zz, th0 * zz * 1e3, color=INK2, ls=sty["ls"], label=f"{E} GeV: $\\theta_0$ = {th0 * 1e3:.2f} mrad")
    ax[0].set_xlabel("distance of the 0.2 X0 from the detector [m]")
    ax[0].set_ylabel("R90 of layers 1-2, sim $-$ sim_v7 [mm]")
    ax[0].set_title("simulated widening (red); blue lines: what the data need", fontsize=10)
    ax[0].legend(fontsize=8)
    ax[1].axhspan(5.5 / 2, 5.5, color=GRID, alpha=0.8, lw=0)
    ax[1].text(1, 4.3, "half a pad .. one pad (5.5 mm)", fontsize=8, color=INK2)
    ax[1].set_xlabel("distance L [m]")
    ax[1].set_ylabel("$\\theta_0 \\cdot L$ [mm]")
    ax[1].set_title("electron-photon separation from multiple scattering in 0.2 X0", fontsize=10)
    ax[1].legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(f"{out}/lever.png", dpi=120)
    plt.close(fig)


def fig_delay(out, test_dir):
    import uproot
    fig, ax = plt.subplots(1, 2, figsize=(11, 4.2))
    z = np.linspace(-70, 0, 300)
    t = (-z * 1e3) / 299.792458                       # ns from origin to the detector front (z ~ 0)
    ax[0].axhspan(0, 200, color=GRID, alpha=0.7, lw=0)        # neutral: blue is reserved for the data
    ax[0].text(-69, 185, "fast-shaper search window (0 .. 200 ns)", fontsize=8, color=INK2, va="top")
    ax[0].plot(z, t, color=INK2, ls="--", label="arrival time, no offset")
    ax[0].plot(z, t - np.maximum(0, -z * 1e3 - 2650) / 299.792458, color=ps.DIGI,
               label="arrival time $-$ TimeOffsetNs (digitiser)")
    for zz, lab in ((-2.65, "sim_v7: -2.65 m, 8.8 ns"), (-65.5, "new: -65.5 m, 218 ns")):
        tt = -zz * 1e3 / 299.792458
        ax[0].plot([zz], [tt], "o", color=INK, ms=6)
        ax[0].annotate(lab, (zz, tt), xytext=(8 if zz < -10 else -150, 8), textcoords="offset points", fontsize=8)
    ax[0].set_xlabel("beam origin z [m]")
    ax[0].set_ylabel("time of the first energy deposit [ns]")
    ax[0].set_ylim(0, 240)
    ax[0].legend(fontsize=8, loc="center left")
    ax[0].set_title("why the H2-like line needs the offset", fontsize=10)
    if test_dir:
        tb = np.linspace(0, 60, 121)
        for (d, lab), sty in zip(((f"{test_dir}/t3/d_0", "sim_v7 geometry, beam at -2.65 m"),
                                  (f"{test_dir}/t2/d_0", "new geometry, beam at -65.5 m, TimeOffsetNs = 209.65 ns")),
                                 ps.reds(2)):
            p = f"{d}/ecal.root"
            if not os.path.exists(p):
                continue
            a = uproot.open(p)["ecal"].arrays(["hit_time", "hit_energy"], library="np")
            ht = np.concatenate(a["hit_time"]); he = np.concatenate(a["hit_energy"])
            ht = ht[(he >= 0.5) & (ht >= 0)]
            nev = len(a["hit_time"])
            h, _ = np.histogram(ht, bins=tb, density=True)
            ps.sim_hist(ax[1], h, tb, sty, label=f"{lab} ({len(ht) / nev:.0f} hits/ev)")
        p = f"{test_dir}/t2/d_nooff/ecal.root"
        if os.path.exists(p):
            a = uproot.open(p)["ecal"].arrays(["hit_energy"], library="np")
            n = sum(int((x >= 0.5).sum()) for x in a["hit_energy"]) / max(len(a["hit_energy"]), 1)
            ax[1].text(0.98, 0.55, f"new geometry WITHOUT the offset:\n{n:.1f} hits/ev", ha="right",
                       transform=ax[1].transAxes, fontsize=8.5, color=INK)
    ax[1].set_xlabel("trigger time of the hit [ns]")
    ax[1].set_ylabel("hits (normalised)")
    ax[1].set_title("74 GeV e$^-$: trigger time with the offset, far vs near", fontsize=10)
    ax[1].legend(fontsize=7.5, loc="upper right")
    fig.tight_layout()
    fig.savefig(f"{out}/delay.png", dpi=120)
    plt.close(fig)


# ------------------------------------------------------------------------------------------------------------- s2
def read_cm(run):
    rows = []
    for line in open(f"{st.OUT}/occupancy/cm/cm_trigger_{run}.txt"):
        m = re.match(r"\s*(\d+)-\s*(\d+)\s+(-?[\d.]+) \|\s+([\d.]+)\s+([\d.]+)\s+([\d.]+) \|\s+([\d.]+)", line)
        if m:
            rows.append([float(x) for x in m.groups()])
    return np.array(rows)


def fig_pedestal(out):
    fig, ax = plt.subplots(1, 2, figsize=(11, 4.2))
    qq = np.logspace(1.5, 4.3, 100)
    # data only: one blue per run, told apart by the marker
    for (run, lab), mk, col in zip(CM_RUNS.items(), ("o", "s", "^"), (ps.DATA, "#1f78a8", "#7cc7ea")):
        r = read_cm(run)
        lo, hi = r[:, 0], np.minimum(r[:, 1], 16000)
        qc = np.where(lo > 0, np.sqrt(np.maximum(lo, 1) * hi), 100)
        qe = [qc - np.maximum(lo, 30), hi - qc]          # the first bin starts at 0: drawn from 30 on the log axis
        ps.data_points(ax[0], qc, r[:, 2], None, None, label=f"data {lab} ({run})", color=col, fmt=mk, xerr=qe)
        ps.data_points(ax[1], qc, r[:, 6], None, None, label=f"data {lab}", color=col, fmt=mk, xerr=qe)
    qq = np.logspace(1.5, np.log10(1.4e4), 100)              # the measured range only
    # the pedestal measurement is a straight line WITH an offset (fit of a + kQ to the three runs together); the
    # digitiser uses k Q through the origin, with k fitted to the hit count of the th230 simulation
    for a0, k, ls, what in ((-0.853, -0.953e-3, ":", "sag measured in pedestals: B = -0.85 - 0.95e-3 $\\cdot$ Q"),
                            (0.0, -1.5e-3, "--", "digitiser, k fitted to the hit count (th230): B = -1.5e-3 $\\cdot$ Q")):
        ax[0].plot(qq, a0 + k * qq, ls, color=ps.DIGI, lw=1.4, label=what)
    ax[0].set_xscale("log")
    ax[0].set_xlabel("Q = HG-ADC charge of the other channels of the chip (same SCA) [ADC]")
    ax[0].set_ylabel("<B>: baseline of the channel [ADC]")
    ax[0].set_title("measured baseline sag in electron runs (channels without signal)", fontsize=10)
    ax[0].legend(fontsize=7.5)
    ax[1].set_xscale("log")
    ax[1].set_xlabel("Q [ADC]")
    ax[1].set_ylabel("trigger turn-on mean, in true amplitude [ADC]")
    ax[1].set_title("the discriminator sees the sagged baseline: the effective threshold rises", fontsize=10)
    ax[1].legend(fontsize=7.5)
    fig.tight_layout()
    fig.savefig(f"{out}/pedestal.png", dpi=120)
    plt.close(fig)


def fig_kscan(out, sets=None):
    """One figure per threshold set, kscan_<set>.png: chip-multiplicity chi2 / chi2(k=0) and hits sim/data against k,
    for that set's groups in occupancy_fit.txt. The sag is fitted per set (a single k does not hold for th210)."""
    lines = open(f"{st.OUT}/occupancy/fit/occupancy_fit.txt").read().splitlines()
    groups = [l for l in lines if l.strip().startswith("model")][0].split()[1:]
    ks, chi, hits = [], [], []
    for l in lines:
        m = re.match(r"lin k=(-?[\d.]+)\s+(.*)", l)
        if m:
            v = np.array(m.group(2).split(), float).reshape(len(groups), 4)
            ks.append(float(m.group(1))); chi.append(v[:, 0]); hits.append(v[:, 1])
    ks, chi, hits = np.array(ks), np.array(chi), np.array(hits)
    by_set = {}
    for j, g in enumerate(groups):
        by_set.setdefault(g.split("_")[0], []).append(j)
    for sname, idx in by_set.items():
        if sets and sname not in sets:
            continue
        idx = sorted(idx, key=lambda j: float(re.search(r"_e([\d.]+)", groups[j]).group(1)))
        fig, ax = plt.subplots(1, 2, figsize=(11, 4.2))
        for j, sty in zip(idx, ps.reds(len(idx))):
            E = re.search(r"_e([\d.]+)", groups[j]).group(1)
            ps.sim_line(ax[0], ks, chi[:, j] / chi[0, j], sty, label=f"sim {E} GeV")
            ps.sim_line(ax[1], ks, hits[:, j], sty, label=f"sim {E} GeV")
        for a in ax:
            a.axvline(-1.5, color=INK, lw=1, ls="--")
            a.axvline(-0.95, color=INK2, lw=1, ls=":")
            a.set_xlabel("k [10$^{-3}$ ADC per ADC]")
        ax[0].text(-1.45, 1.6, "chosen\n-1.5e-3", fontsize=8)
        ax[0].text(-0.9, 1.6, "pedestals\n-0.95e-3", fontsize=8, color=INK2)
        ax[0].set_ylabel("chip-multiplicity $\\chi^2$ / $\\chi^2$(k = 0)")
        ax[0].set_ylim(0, 2)
        ax[0].set_title(f"{sname}: shape of the hits-per-chip spectrum", fontsize=10)
        ax[1].axhline(1, color=INK2, lw=1)
        ax[1].set_ylabel("hits / event, sim / data")
        ax[1].set_title(f"{sname}: total hit count (H2-like geometry)", fontsize=10)
        ax[1].legend(fontsize=7.5)
        fig.tight_layout()
        fig.savefig(f"{out}/kscan_{sname}.png", dpi=120)
        plt.close(fig)


def load_ecal(path):
    import uproot
    br = {"e": "hit_energy", "hg": "hit_hg", "slab": "hit_slab", "chip": "hit_chip", "chan": "hit_chan",
          "x": "hit_x", "y": "hit_y", "lg": "hit_lg"}
    a = uproot.open(path)["ecal"].arrays(list(br.values()), library="np")
    n = len(a["hit_energy"])
    s = {"ev": np.concatenate([np.full(len(x), i) for i, x in enumerate(a["hit_energy"])]), "nev": n}
    for k, b in br.items():
        s[k] = np.concatenate(a[b])
    s["ix"] = s["iy"] = np.zeros(len(s["e"]), int)
    s["slab"] = s["slab"].astype(int); s["chip"] = s["chip"].astype(int); s["e"] = s["e"].astype(np.float32)
    return st.subset(s, s["e"] >= 0.5)


def fig_impl(out, test_dir):
    s0 = load_ecal(f"{test_dir}/t2/d_0/ecal.root")
    s1 = load_ecal(f"{test_dir}/t2/d_on/ecal.root")
    ph = [st.apply_baseline(s0, "th230", -1.5e-3, np.random.default_rng(i)) for i in range(5)]
    edges = np.array([0, 250, 500, 1000, 2000, 4000, 8000, 30000.0])
    q0 = np.digitize(st.chip_charge(s0, "adc"), edges) - 1

    def per_bin(s):
        # chip charge of the k=0 sample, per (event, slab, chip): bins the SAME chips in every variant
        key0 = st.chip_key(s0)
        qmap = dict(zip(key0, q0))
        b = np.array([qmap.get(k, -1) for k in st.chip_key(s)])
        return np.array([(b == i).sum() for i in range(len(edges) - 1)]), \
            np.array([s["e"][b == i].sum() for i in range(len(edges) - 1)])

    n0, e0 = per_bin(s0)
    n1, e1 = per_bin(s1)
    nph = np.array([per_bin(p)[0] for p in ph])
    eph = np.array([per_bin(p)[1] for p in ph])
    qc = np.sqrt(np.maximum(edges[:-1], 60) * edges[1:])
    fig, ax = plt.subplots(1, 2, figsize=(11, 4.2))
    ok = n0 > 50
    s_dig, s_ph = ps.reds(2)[::-1]
    ps.sim_line(ax[0], qc[ok], (n1 / n0)[ok], s_dig, label="AdcDigitizer OccupancySlope = -1.5e-3")
    ax[0].errorbar(qc[ok], (nph / n0).mean(0)[ok], (nph / n0).std(0)[ok], fmt="s", ls=s_ph["ls"], color=s_ph["color"],
                   ms=4, lw=1.4, label="post-hoc sim_tuning.apply_baseline (5 seeds)")
    ax[0].set_xscale("log")
    ax[0].set_xlabel("Q of the hit's chip, before the sag [ADC]")
    ax[0].set_ylabel("hits kept / hits at k = 0")
    ax[0].set_title("74 GeV e$^-$, 200 events, same Geant4 events", fontsize=10)
    ax[0].legend(fontsize=8)
    r1 = (e1 / n1) / (e0 / n0)
    rph = (eph / nph) / (e0 / n0)
    ps.sim_line(ax[1], qc[ok], r1[ok], s_dig, label="AdcDigitizer")
    ax[1].errorbar(qc[ok], rph.mean(0)[ok], rph.std(0)[ok], fmt="s", ls=s_ph["ls"], color=s_ph["color"], ms=4, lw=1.4,
                   label="post-hoc")
    ax[1].set_xscale("log")
    ax[1].set_xlabel("Q [ADC]")
    ax[1].set_ylabel("E/hit / E/hit at k = 0")
    ax[1].set_title("the lost hits are the small ones: E/hit goes UP", fontsize=10)
    ax[1].legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(f"{out}/impl_check.png", dpi=120)
    plt.close(fig)
    tot = (len(s1["e"]) / len(s0["e"]), np.mean([len(p["e"]) for p in ph]) / len(s0["e"]),
           s1["e"].mean() / s0["e"].mean(), np.mean([p["e"].mean() for p in ph]) / s0["e"].mean())
    txt = (f"74 GeV, {int(s0['nev'])} events: hits ratio digitiser {tot[0]:.4f}, post-hoc {tot[1]:.4f}; "
           f"E/hit ratio digitiser {tot[2]:.4f}, post-hoc {tot[3]:.4f}")
    open(f"{out}/impl_check.txt", "w").write(txt + "\n")
    print(txt)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("outdir")
    ap.add_argument("--test-dir", default=None)
    ap.add_argument("--only", default="s1,s2")
    args = ap.parse_args()
    os.makedirs(args.outdir, exist_ok=True)
    only = args.only.split(",")
    if "s1" in only:
        fig_profile_long(args.outdir); fig_lateral(args.outdir); fig_lever(args.outdir)
        fig_delay(args.outdir, args.test_dir)
    if "s2" in only:
        fig_pedestal(args.outdir); fig_kscan(args.outdir)
        if args.test_dir:
            fig_impl(args.outdir, args.test_dir)


if __name__ == "__main__":
    main()

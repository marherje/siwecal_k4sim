#!/usr/bin/env python
"""
Common tools of the sim_tuning study (2026-10-16): hit-level samples of data and simulation, the observables every
strategy is judged on, detector heat maps, and the offline emulations of the candidate digitiser changes.

A sample is the shower-like events (nhit > 0.5 x p90, as final_panels.py) of one or more valtrees (0.5 MIP hit cut of
the PID stage, identical in data and simulation), flattened to one row per hit:
    ev slab chip chan ix iy x y e hg lg      + per event: bar_x bar_y
(ix, iy) = cell column/row on the 32 x 32 pad grid, from hit_x/hit_y.

Observables (observables()): hits/event and energy/event (Gaussian core), E/hit (mean over hits), per layer hits and
E/hit, per chip (slab x 16) and per cell (slab x 32 x 32) hits/event and energy/event, radial profile about the
event barycentre. compare() turns two of them into sim/data ratios and the scalar figures of merit:
    hits, E, E/hit            global sim/data
    layer_chi2                per-layer hits/event, (s-d)^2/(var_s+var_d) per layer, mean
    chip_rms, cell_rms        hit-weighted RMS of log2(sim/data) over chips / cells with >= MIN_HITS hits in both
                              (cells: after subtracting the 5x5 median, i.e. the local residual)

Emulations (apply_*) act on a simulation sample and return a new one; they never touch the data.
"""
import glob
import os
import sys

import numpy as np
import uproot
import yaml
from scipy.special import erf
from scipy.ndimage import median_filter

_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(_HERE))
from analysis.compare_event_level import gauss_core  # noqa: E402
from analysis.final_panels import load_groups, RECO  # noqa: E402

B = "/eos/experiment/drdcalo/siw-ecal/TB2026-06"
V7 = f"{B}/Simulation/Processed/final_v7"
OUT = f"{B}/Simulation/tests/sim_tuning"
NL, NCHIP, NG = 15, 16, 32
PITCH = 172.64 / 31.0
X0 = -86.32
MIN_HITS = 25
ADC_PER_MIP = 18.75
LG_SWITCH = 1500.0
BRANCHES = ["nhit_chan", "bar_x", "bar_y", "hit_slab", "hit_chip", "hit_chan", "hit_x", "hit_y", "hit_energy",
            "hit_hg", "hit_lg"]
R_EDGES = np.array([0, 5.5, 11, 16.5, 22, 33, 44, 66, 99, 200.0])


# ---------------------------------------------------------------------------------------------------- samples
def data_paths(group):
    return [p for r in group["runs"] for p in sorted(glob.glob(f"{RECO}/{r}/*.valtree.root"))]


def sim_paths(name, campaign=V7, kind="digi"):
    return sorted(glob.glob(f"{campaign}/pid/{name}_{kind}/*.valtree.root"))


def load(paths, nmax=40000, cache=None, alone=False):
    """Shower-like events of the valtrees, flattened per hit. Cached as .npz when cache is given.
    alone=True (data): keep only events alone in their DAQ acquisition (single_run of the reconstruction's ecal tree,
    <run>.root next to <run>.valtree.root, matched entry by entry); the shower-like cut is set on all events."""
    if cache and os.path.exists(cache):
        return dict(np.load(cache))
    parts, single, left = [], [], nmax
    for p in paths:
        if left <= 0:
            break
        a = uproot.open(p)["ecal"].arrays(BRANCHES + (["event", "spill"] if alone else []), entry_stop=left,
                                          library="np")
        if alone:
            r = uproot.open(p.replace(".valtree.root", ".root"))["ecal"].arrays(
                ["event", "spill", "single_run"], entry_stop=left, library="np")
            if not (np.array_equal(r["event"], a["event"]) and np.array_equal(r["spill"], a["spill"])):
                raise RuntimeError(f"{p}: valtree and ecal tree entries do not match")
            single.append(r["single_run"].astype(bool))
        left -= len(a["nhit_chan"])
        parts.append(a)
    a = {k: np.concatenate([p[k] for p in parts]) for k in BRANCHES}
    keep = a["nhit_chan"] > 0.5 * np.percentile(a["nhit_chan"], 90)
    if alone:
        keep &= np.concatenate(single)
    sel = np.nonzero(keep)[0]
    cnt = np.array([len(a["hit_energy"][i]) for i in sel])
    s = {"ev": np.repeat(np.arange(len(sel), dtype=np.int32), cnt), "nev": np.array(len(sel)),
         "bar_x": a["bar_x"][sel].astype(np.float32), "bar_y": a["bar_y"][sel].astype(np.float32)}
    for k, b in (("slab", "hit_slab"), ("chip", "hit_chip"), ("chan", "hit_chan"), ("x", "hit_x"), ("y", "hit_y"),
                 ("e", "hit_energy"), ("hg", "hit_hg"), ("lg", "hit_lg")):
        s[k] = np.concatenate([a[b][i] for i in sel]) if len(sel) else np.zeros(0)
    for k in ("slab", "chip", "chan"):
        s[k] = s[k].astype(np.int16)
    for k in ("x", "y", "e", "hg", "lg"):
        s[k] = s[k].astype(np.float32)
    s["ix"] = np.clip(np.rint((s["x"] - X0) / PITCH), 0, NG - 1).astype(np.int16)
    s["iy"] = np.clip(np.rint((s["y"] - X0) / PITCH), 0, NG - 1).astype(np.int16)
    if cache:
        os.makedirs(os.path.dirname(cache), exist_ok=True)
        np.savez(cache, **s)
    return s


def subset(s, keep, e=None):
    """New sample with the hits in mask keep (and new energies e for all hits, before masking)."""
    out = {k: v for k, v in s.items() if k in ("nev", "bar_x", "bar_y")}
    for k in ("ev", "slab", "chip", "chan", "x", "y", "hg", "lg", "ix", "iy"):
        out[k] = s[k][keep]
    out["e"] = (s["e"] if e is None else e)[keep]
    return out


def chip_key(s):
    return (s["ev"].astype(np.int64) * NL + s["slab"]) * NCHIP + s["chip"]


def chip_charge(s, unit="mip"):
    """Q of each hit: summed energy of the OTHER hits of its chip in the event (MIP, or HG ADC with unit='adc')."""
    _, inv = np.unique(chip_key(s), return_inverse=True)
    q = np.bincount(inv, weights=s["e"])[inv] - s["e"]
    return q * ADC_PER_MIP if unit == "adc" else q


def chip_nhits(s):
    _, inv, cnt = np.unique(chip_key(s), return_inverse=True, return_counts=True)
    return cnt[inv]


# ---------------------------------------------------------------------------------------------------- observables
def observables(s):
    nev = int(s["nev"])
    ev, e, sl = s["ev"], s["e"].astype(float), s["slab"]
    nh_ev = np.bincount(ev, minlength=nev).astype(float)
    E_ev = np.bincount(ev, weights=e, minlength=nev)
    o = {"nev": nev, "hits": gauss_core(nh_ev)[0], "E": gauss_core(E_ev)[0], "eph": e.sum() / max(len(e), 1),
         "hits_mean": nh_ev.mean(), "E_mean": E_ev.mean()}
    # per layer: mean and variance of the per-event count
    lay = np.zeros((nev, NL)); np.add.at(lay, (ev, sl), 1)
    layE = np.zeros((nev, NL)); np.add.at(layE, (ev, sl), e)
    o["layer_hits"], o["layer_hits_var"] = lay.mean(0), lay.var(0) / nev
    o["layer_E"] = layE.mean(0)
    o["layer_eph"] = layE.sum(0) / np.maximum(lay.sum(0), 1)
    ci = sl.astype(int) * NCHIP + s["chip"]
    o["chip_hits"] = np.bincount(ci, minlength=NL * NCHIP).reshape(NL, NCHIP) / nev
    o["chip_E"] = np.bincount(ci, weights=e, minlength=NL * NCHIP).reshape(NL, NCHIP) / nev
    cc = (sl.astype(int) * NG + s["iy"]) * NG + s["ix"]
    o["cell_hits"] = np.bincount(cc, minlength=NL * NG * NG).reshape(NL, NG, NG) / nev
    o["cell_E"] = np.bincount(cc, weights=e, minlength=NL * NG * NG).reshape(NL, NG, NG) / nev
    r = np.hypot(s["x"] - s["bar_x"][ev], s["y"] - s["bar_y"][ev])
    rb = np.clip(np.digitize(r, R_EDGES) - 1, 0, len(R_EDGES) - 2)
    o["radial_E"] = np.bincount(sl.astype(int) * (len(R_EDGES) - 1) + rb, weights=e,
                                minlength=NL * (len(R_EDGES) - 1)).reshape(NL, -1) / nev
    o["radial_hits"] = np.bincount(sl.astype(int) * (len(R_EDGES) - 1) + rb,
                                   minlength=NL * (len(R_EDGES) - 1)).reshape(NL, -1) / nev
    return o


def log2ratio(sim, dat, ns, nd, min_hits=MIN_HITS):
    """log2(sim/data) of per-event maps, NaN where either sample has < min_hits hits in total."""
    ok = (sim * ns >= min_hits) & (dat * nd >= min_hits)
    with np.errstate(divide="ignore", invalid="ignore"):
        return np.where(ok, np.log2(sim / dat), np.nan)


def local_residual(m):
    """m minus the median of its 5x5 neighbourhood, per slab (NaN-aware: NaN cells are ignored)."""
    out = np.full_like(m, np.nan)
    for l in range(m.shape[0]):
        a = m[l]
        fill = np.where(np.isnan(a), np.nanmedian(a) if np.isfinite(a).any() else 0, a)
        out[l] = a - median_filter(fill, size=5, mode="nearest")
    return out


def wrms(x, w):
    k = np.isfinite(x) & (w > 0)
    return float(np.sqrt(np.average(x[k] ** 2, weights=w[k]))) if k.any() else np.nan


def compare(so, do):
    """sim/data figures of merit of two observables() dicts."""
    ns, nd = so["nev"], do["nev"]
    c = {"hits": so["hits"] / do["hits"], "E": so["E"] / do["E"], "eph": so["eph"] / do["eph"]}
    c["layer_ratio"] = so["layer_hits"] / np.maximum(do["layer_hits"], 1e-9)
    c["layer_eph_ratio"] = so["layer_eph"] / np.maximum(do["layer_eph"], 1e-9)
    v = so["layer_hits_var"] + do["layer_hits_var"]
    c["layer_chi2"] = float(np.mean((so["layer_hits"] - do["layer_hits"]) ** 2 / np.maximum(v, 1e-12)))
    c["chip_log2"] = log2ratio(so["chip_hits"], do["chip_hits"], ns, nd)
    eph_s = so["chip_E"] / np.maximum(so["chip_hits"], 1e-12)
    eph_d = do["chip_E"] / np.maximum(do["chip_hits"], 1e-12)
    with np.errstate(divide="ignore", invalid="ignore"):
        c["chip_eph_log2"] = np.where(np.isfinite(c["chip_log2"]), np.log2(eph_s / eph_d), np.nan)
    c["cell_log2"] = log2ratio(so["cell_hits"], do["cell_hits"], ns, nd)
    c["cell_res"] = local_residual(c["cell_log2"])
    c["chip_rms"] = wrms(c["chip_log2"], do["chip_hits"])
    c["chip_eph_rms"] = wrms(c["chip_eph_log2"], do["chip_hits"])
    c["cell_rms"] = wrms(c["cell_log2"], do["cell_hits"])
    c["cell_res_rms"] = wrms(c["cell_res"], do["cell_hits"])
    c["radial_ratio"] = so["radial_E"] / np.maximum(do["radial_E"], 1e-12)
    return c


# ---------------------------------------------------------------------------------------------------- calibration
def trigger_params(th, path=os.path.join(os.path.dirname(_HERE), "mappings", "digi_calibration.yml")):
    """(mu, sigma) of the per-slab trigger turn-on in HG ADC, arrays over the 15 slabs."""
    c = yaml.safe_load(open(path))["thresholds"][th]
    mu = np.full(NL, c["threshold_adc"]); sg = np.full(NL, c["sigma_mip"] * c["adc_per_mip"])
    for l, o in (c.get("slab_overrides") or {}).items():
        mu[int(l)] = o["threshold_adc"]; sg[int(l)] = o.get("sigma_adc", sg[int(l)])
    return mu, sg


def p_trig(a, mu, sg):
    return 0.5 * (1 + erf((a - mu) / (np.sqrt(2) * sg)))


# ---------------------------------------------------------------------------------------------------- emulations
def apply_baseline(s, th, slope, rng, slope_map=None, q0_adc=0.0):
    """2.1: chip baseline drop B = slope x (Q - Q0) [HG ADC], Q = HG-ADC charge of the other hits of the chip.
    HG-read hits lose B in amplitude and energy; every hit survives the trigger with P(a + B) / P(a) (per-slab erf
    turn-on) and the 0.5 MIP cut. slope_map[slab, chip] (ADC per ADC) replaces the scalar where finite."""
    q = chip_charge(s, "adc")
    k = np.full(len(q), slope)
    if slope_map is not None:
        m = slope_map[s["slab"], s["chip"]]
        k = np.where(np.isfinite(m), m, k)
    b = k * np.maximum(q - q0_adc, 0.0)                                        # <= 0
    mu, sg = trigger_params(th)
    mu, sg = mu[s["slab"]], sg[s["slab"]]
    hg = s["hg"].astype(float)
    surv = np.clip(p_trig(hg + b, mu, sg) / np.maximum(p_trig(hg, mu, sg), 1e-9), 0, 1)
    hgread = hg < LG_SWITCH
    e = np.where(hgread, s["e"] + b / ADC_PER_MIP, s["e"])
    keep = (rng.random(len(e)) < surv) & (e >= 0.5)
    return subset(s, keep, e.astype(np.float32))


def apply_lowgain(s, q_c, r_q, r_ref, switch=LG_SWITCH):
    """2.3: hits read in low gain (HG >= switch) scaled by r(Q)/r_ref, r = measured LG/HG ratio against chip charge
    Q [MIP] (bin centres q_c, values r_q)."""
    q = chip_charge(s)
    r = np.interp(np.log(np.maximum(q, 1.0)), np.log(q_c), r_q)
    e = np.where(s["hg"] >= switch, s["e"] * r / r_ref, s["e"]).astype(np.float32)
    return subset(s, np.ones(len(e), bool), e)


def apply_efficiency_map(s, eff, rng):
    """2.2: per-cell (or per-chip, broadcast) efficiency eff[slab, iy, ix] <= 1: drop sim hits with prob 1 - eff."""
    p = np.clip(eff[s["slab"], s["iy"], s["ix"]], 0, 1)
    p = np.where(np.isfinite(p), p, 1.0)
    return subset(s, rng.random(len(p)) < p)


def apply_chip_dropout(s, p_drop, rng):
    """2.2: whole chips missing from an event with probability p_drop[slab, chip] (readout dropouts)."""
    key = chip_key(s)
    u, inv = np.unique(key, return_inverse=True)
    sl, ch = (u // NCHIP) % NL, u % NCHIP
    gone = rng.random(len(u)) < np.nan_to_num(p_drop[sl, ch])
    return subset(s, ~gone[inv])

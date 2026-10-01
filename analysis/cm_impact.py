#!/usr/bin/env python
"""
Impact study, before implementing it, of the chip common-mode baseline drop (the median ADC - pedestal of the channels
WITHOUT a hit bit falls with the number of hit bits in the SCA: -0.8 ADC at 0-3 bits to -5.3 ADC at 20-27, measured
on run 286 raw data with the fixed pedestals).

Case 1 (correct the data): every data hit gets +(sag(N) - sag(N0)) ADC back, N = hits of its chip in the event,
        N0 the muon-calibration occupancy (the pedestal and MIP tables already contain sag(N0)).
Case 2 (model it in the simulation): every simulated hit loses (sag(N) - sag(N0)) ADC; it is lost if it no longer
        passes the trigger (erf turn-on, threshold/width of the set) -- drawn with the survival probability
        P(a')/P(a) -- or the 0.5 MIP cut.

Both on the valtrees; per group: nhit, energy [MIP] (Gaussian core), E/hit (median) and halo hits (15-60 mm,
layers 3-8), as sim / data for: as now, case 1, case 2, both.

    python analysis/cm_impact.py <out.txt> [group names]
"""
import sys

import numpy as np
import uproot

sys.path.insert(0, __file__.rsplit("/analysis/", 1)[0])
sys.path.insert(0, __file__.rsplit("/siwecal_k4sim/", 1)[0] + "/siwecal-tb2026")
from analysis.compare_event_level import gauss_core  # noqa: E402
from analysis.final_panels import load_groups, sim_valtree, RECO  # noqa: E402
from scipy.special import erf  # noqa: E402

# measured baseline (ADC - pedestal) of the no-bit channels against the hit bits in the SCA (run 286, fixed data)
SAG_N = np.array([1.5, 5.5, 9.5, 13.5, 17.5, 21.5, 25.5, 29.5, 33.5, 37.5, 41.5, 45.5, 49.5])
SAG_V = np.array([-0.75, -1.82, -3.18, -4.36, -5.04, -5.33, -5.34, -5.10, -4.44, -3.98, -2.98, 1.11, 3.97])
N0 = 1.5
TRIG = {"th210": (16.11, 3.33), "th220": (20.23, 3.58), "th230": (26.07, 3.74)}
ADC_PER_MIP = 18.75


def sag(n):
    return np.interp(n, SAG_N, SAG_V) - np.interp(N0, SAG_N, SAG_V)   # <= 0: ADC lost relative to calibration


def p_trig(a, th):
    mu, sg = TRIG[th]
    return 0.5 * (1 + erf((a - mu) / (np.sqrt(2) * sg)))


def load(path, nmax=20000):
    a = uproot.open(path)["ecal"].arrays(["nhit_chan", "bar_x", "bar_y", "hit_x", "hit_y", "hit_slab", "hit_chip",
                                          "hit_hg", "hit_energy"], entry_stop=nmax, library="np")
    k = np.nonzero(a["nhit_chan"] > 0.5 * np.percentile(a["nhit_chan"], 90))[0][:8000]
    return {b: a[b][k] for b in a}


def observables(ev, mode, th, rng):
    """mode: 'none', 'restore' (case 1, data) or 'sag' (case 2, sim)."""
    nh, E, eph, halo = [], [], [], []
    for i in range(len(ev["nhit_chan"])):
        s, c, hg, e = ev["hit_slab"][i], ev["hit_chip"][i], ev["hit_hg"][i].astype(float), ev["hit_energy"][i].astype(float)
        if mode != "none":
            key = s.astype(int) * 16 + c.astype(int)
            _, inv, cnt = np.unique(key, return_inverse=True, return_counts=True)
            d = sag(cnt[inv])                                   # <= 0
            mip_per_adc = np.where(hg > 20, e / np.maximum(hg, 1e-9), 1.0 / ADC_PER_MIP)
            if mode == "restore":
                e = e - d * mip_per_adc
            else:
                new = hg + d
                keep = rng.random(len(hg)) < np.clip(p_trig(new, th) / np.maximum(p_trig(hg, th), 1e-9), 0, 1)
                e = e + d * mip_per_adc
                keep &= e >= 0.5
                s, e = s[keep], e[keep]
                hx, hy = ev["hit_x"][i][keep], ev["hit_y"][i][keep]
        if mode != "sag":
            hx, hy = ev["hit_x"][i], ev["hit_y"][i]
        r = np.hypot(hx - ev["bar_x"][i], hy - ev["bar_y"][i])
        nh.append(len(e)); E.append(e.sum()); eph.append(e.sum() / max(len(e), 1))
        halo.append(((r > 15) & (r < 60) & (s >= 3) & (s <= 8)).sum())
    nh, E = np.array(nh, float), np.array(E)
    return {"nhit": gauss_core(nh)[0], "E": gauss_core(E)[0], "E/hit": np.median(eph), "halo": np.mean(halo)}


def main():
    out = sys.argv[1]; only = set(sys.argv[2:])
    rng = np.random.default_rng(1)
    lines = [f"{'group':<22}{'case':<12}" + "".join(f"{k:>9}" for k in ("nhit", "E", "E/hit", "halo")) + "   (sim / data)"]
    for g in load_groups():
        if only and g["name"] not in only:
            continue
        th = g["set"]
        d = load(f"{RECO}/{g['runs'][0]}/ecal_{g['runs'][0]}.valtree.root")
        s = load(sim_valtree(g, "digi"))
        D0, D1 = observables(d, "none", th, rng), observables(d, "restore", th, rng)
        S0, S2 = observables(s, "none", th, rng), observables(s, "sag", th, rng)
        for lab, sv, dv in (("as now", S0, D0), ("case 1", S0, D1), ("case 2", S2, D0), ("both", S2, D1)):
            lines.append(f"{g['name']:<22}{lab:<12}" + "".join(f"{sv[k] / dv[k]:9.3f}" for k in ("nhit", "E", "E/hit", "halo")))
        print("\n".join(lines[-4:]), flush=True)
    open(out, "w").write("\n".join(lines) + "\n")


if __name__ == "__main__":
    main()

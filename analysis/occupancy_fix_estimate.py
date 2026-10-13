#!/usr/bin/env python
"""
How much of the sim/data difference would the two measured occupancy effects remove? Estimate on the simulated hits,
before implementing them in the digitiser.

Both are functions of the charge Q of the other hits of the same chip (one trigger per chip and event in the sim):
  baseline   B(Q) = -0.80 ADC per 1000 ADC of Q (cm_trigger_test.py / cm_per_chip.py, retrigger veto applied), relative
             to the muon calibration. Every high-gain hit loses B; it survives the trigger with probability
             P(a + B) / P(a) (erf turn-on of the threshold set, as in cm_impact.py) and the 0.5 MIP cut.
  low gain   hits read in low gain (high gain - pedestal > 1500 ADC, the event builder's switch) are scaled by
             r(Q) / r_ref, r(Q) = the measured LG/HG ratio against Q (cm_mip_peak.py, 74 GeV). The LG->HG anchor was
             fitted on showers (eudaq 166, 99 GeV), so it already holds an average Q: r_ref is r averaged over the
             hits in the anchor's fit range (150-1500 ADC) of the highest-energy sample. Upper bound: r_ref = quiet chip.

    python analysis/occupancy_fix_estimate.py <out.txt> [group names]

Writes sim/data of hits, energy (Gaussian core) and E/hit (median) for: as now, + baseline, + low gain, both, and both
with the quiet-chip reference.
"""
import sys

import numpy as np

from cm_impact import load, p_trig, ADC_PER_MIP, gauss_core, load_groups, sim_valtree, RECO

B_SLOPE = -0.80e-3
LG_SWITCH = 1500.0
Q_C = np.array([2.2, 8.7, 21, 42, 85, 173, 354, 700])                                # MIP, bin centres
R_Q = np.array([0.1116, 0.1114, 0.1099, 0.1080, 0.1056, 0.1038, 0.1017, 0.1000])     # LG/HG, 74 GeV, LG baseline-corrected
GROUPS = ["th210_S1_e20", "th210_S1_e74", "th230_E4_e7.5", "th230_E4_e10", "th230_E4_e20", "th230_E4_e34",
          "th230_E4_e52", "th230_E4_e74"]


def r_of_q(q_mip):
    return np.interp(np.log(np.maximum(q_mip, 1.0)), np.log(Q_C), R_Q)


def chip_q(ev, i):
    """Q of each hit: energy of the other hits of its chip, in MIP."""
    key = ev["hit_slab"][i].astype(int) * 16 + ev["hit_chip"][i].astype(int)
    e = ev["hit_energy"][i].astype(float)
    _, inv = np.unique(key, return_inverse=True)
    return np.bincount(inv, weights=e)[inv] - e


def r_reference(ev):
    rs = []
    for i in range(len(ev["nhit_chan"])):
        hg = ev["hit_hg"][i].astype(float)
        sel = (hg > 150) & (hg < LG_SWITCH)
        rs.append(r_of_q(chip_q(ev, i)[sel]))
    return float(np.mean(np.concatenate(rs)))


def observables(ev, th, rng, baseline=False, lowgain=False, r_ref=None):
    nh, E, eph = [], [], []
    for i in range(len(ev["nhit_chan"])):
        hg, e = ev["hit_hg"][i].astype(float), ev["hit_energy"][i].astype(float)
        if baseline or lowgain:
            q = chip_q(ev, i)
            lg = hg > LG_SWITCH
            keep = np.ones(len(e), bool)
            if baseline:
                b = B_SLOPE * q * ADC_PER_MIP
                mip_per_adc = np.where(hg > 20, e / np.maximum(hg, 1e-9), 1.0 / ADC_PER_MIP)
                keep = rng.random(len(e)) < np.clip(p_trig(hg + b, th) / np.maximum(p_trig(hg, th), 1e-9), 0, 1)
                e = np.where(lg, e, e + b * mip_per_adc)
                keep &= e >= 0.5
            if lowgain:
                e = np.where(lg, e * r_of_q(q) / r_ref, e)
            e = e[keep]
        nh.append(len(e)); E.append(e.sum()); eph.append(e.sum() / max(len(e), 1))
    return {"nhit": gauss_core(np.array(nh, float))[0], "E": gauss_core(np.array(E))[0], "E/hit": np.median(eph)}


def main():
    out = sys.argv[1]
    want = sys.argv[2:] or GROUPS
    rng = np.random.default_rng(1)
    groups = {g["name"]: g for g in load_groups()}
    sims = {n: load(sim_valtree(groups[n], "digi")) for n in want}
    top = max(want, key=lambda n: float(groups[n]["E"]))
    r_ref = r_reference(sims[top])
    lines = [f"LG/HG reference r_ref = {r_ref:.4f} (hits 150-1500 ADC of {top}); quiet chip r = {R_Q[0]:.4f}",
             f"{'group':<16}{'case':<26}" + "".join(f"{k:>9}" for k in ("nhit", "E", "E/hit")) + "   (sim / data)"]
    print(lines[0], flush=True)
    for n in want:
        g = groups[n]
        th = g["set"]
        d = observables(load(f"{RECO}/{g['runs'][0]}/ecal_{g['runs'][0]}.valtree.root"), th, rng)
        s = sims[n]
        cases = [("as now", {}), ("+ baseline", {"baseline": True}), ("+ low gain", {"lowgain": True, "r_ref": r_ref}),
                 ("both", {"baseline": True, "lowgain": True, "r_ref": r_ref}),
                 ("both, LG vs quiet chip", {"baseline": True, "lowgain": True, "r_ref": R_Q[0]})]
        for lab, kw in cases:
            sv = observables(s, th, rng, **kw)
            lines.append(f"{n:<16}{lab:<26}" + "".join(f"{sv[k] / d[k]:9.3f}" for k in ("nhit", "E", "E/hit")))
        print("\n".join(lines[-5:]), flush=True)
    open(out, "w").write("\n".join(lines) + "\n")


if __name__ == "__main__":
    main()

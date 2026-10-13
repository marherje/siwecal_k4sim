#!/usr/bin/env python
"""
Is the gain of a busy chip lower? Position of the 1-MIP peak (high gain) and the low/high gain ratio against the charge
in the chip-SCA, from the decoded data alone.

The baseline drop of a busy chip-SCA (cm_trigger_test.py) is an offset: about -1 ADC per 1000 ADC of charge, i.e. 1-3 %
of the energy per hit. The energy per hit of busy chips is 16-42 % below the simulation, which would need a gain change.
If the amplification of the chip falls with its charge, the 1-MIP peak of the shower halo hits moves down with it.

Acquisitions with fewer than 5 chips (retrigger storms) are skipped. For every channel of every valid chip-SCA
(retrigger columns vetoed, as in cm_trigger_test.py), with or without hit bit,
so the trigger turn-on (which moves with the baseline) does not sculpt the peak:
    s     = (adc_high - pedestal) / MIP_channel                          raw
    s_cor = (adc_high - pedestal - B) / MIP_channel                      B = baseline of the chip-SCA, measured on the
            far non-fired channels when possible, else B0 + b Q (b = -0.80 ADC / 1000 ADC)
binned in Q = charge of the other fired channels of the chip-SCA, in MIP (low gain / k / MIP of each fired channel).
For fired channels with 150 < HG - ped - B < 800 ADC (high gain linear), also the ratio (LG - ped - B_LG) / (HG - ped - B),
with B_LG the low-gain baseline of the same far channels (and without it, ratio_raw), on the chip-SCAs where it is measured.

    python analysis/cm_mip_peak.py <outdir> --run e291 --chunks 10-17 [--max-acq 9000]

Writes cm_mip_peak_<run>.npz (histograms per Q bin). Figures: analysis/cm_mip_peak_plot.py, against the simulation.
"""
import argparse
import os
import sys

import numpy as np
import uproot

sys.path.insert(0, __file__.rsplit("/siwecal_k4sim/", 1)[0] + "/siwecal-tb2026/diagnostics")
from cm_baseline_collect import load_ped, adjacency, run_name, CONV, CAL, MAP, MAP_COB, K_ANCHOR  # noqa: E402

MIPS = CAL.replace("/pedestals/", "/mips/") + "/MIP_pedestalsubmode1_TB2026CERN_run_000th210_highgain.txt"
Q_EDGES = np.array([0, 5, 15, 30, 60, 120, 250, 500, 1e9])        # MIP
S_BINS = np.arange(-1.0, 4.0001, 0.02)                           # MIP
R_BINS = np.arange(0.05, 0.1500001, 0.0005)                      # LG / HG
B_SLOPE, B_QUIET = -0.80e-3, -0.6


def load_mip(path):
    """Per-channel high-gain MIP MPV [ADC], chip median where the channel has no fit."""
    mip = np.full((15, 16, 64), np.nan, np.float32)
    for line in open(path):
        if line.startswith("#"):
            continue
        v = line.split()
        if len(v) >= 4 and float(v[3]) > 0:
            mip[int(v[0]), int(v[1]), int(v[2])] = float(v[3])
    med = np.nanmedian(mip, axis=-1, keepdims=True)
    return np.where(np.isfinite(mip), mip, med)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("outdir")
    ap.add_argument("--run", required=True)
    ap.add_argument("--chunks", default="10-17")
    ap.add_argument("--max-acq", type=int, default=0)
    a = ap.parse_args()
    ped_hg = load_ped(f"{CAL}/Pedestal_TB2026CERN_run_000th210_highgain.txt")
    ped_lg = load_ped(f"{CAL}/Pedestal_TB2026CERN_run_000th210_lowgain.txt")
    mip = load_mip(MIPS)
    _, A2 = adjacency(MAP)
    _, A2c = adjacency(MAP_COB)
    c0, c1 = (int(x) for x in a.chunks.split("-"))
    files = [f"{CONV}/{run_name(a.run)}/chunks/chunk_{i:04d}.root" for i in range(c0, c1 + 1)]
    files = [f for f in files if os.path.exists(f)]
    nq, ns, nr = len(Q_EDGES) - 1, len(S_BINS) - 1, len(R_BINS) - 1
    h = {k: np.zeros(nq * ns, np.int64) for k in ("raw_all", "cor_all", "cor_bit", "cor_nobit")}
    hr = {k: np.zeros(nq * nr, np.int64) for k in ("ratio", "ratio_raw")}
    blg = np.zeros((nq, 2))
    bmeas, n_acq = np.zeros(2), 0
    slabs = np.arange(15)
    for path in files:
        for d in uproot.open(path)["siwecaldecoded"].iterate(
                ["chipid", "nColumns", "bcid", "adc_high", "adc_low", "hitbit_high"], step_size=20, library="np"):
            if a.max_acq and n_acq >= a.max_acq:
                break
            # acquisitions with fewer than 5 chips are retrigger storms (one chip, 15 full columns), not showers
            keep = (d["chipid"] >= 0).sum((1, 2)) >= 5
            if not keep.any():
                continue
            d = {k: v[keep] for k, v in d.items()}
            cid = d["chipid"]
            cidc = np.clip(cid, 0, 15)
            L = slabs[None, :, None]
            P, PL = ped_hg[L, cidc], ped_lg[L, cidc]                             # [a, slab, chip, sca, ch]
            M = mip[L, cidc][..., None, :]                                       # [a, slab, chip, 1, ch]
            adc = d["adc_high"].astype(np.float32)
            hb = d["hitbit_high"].astype(np.float32)
            valid = (np.arange(15)[None, None, None, :] < d["nColumns"][..., None]) & (cid >= 0)[..., None]
            bc = np.where(valid, d["bcid"], -999).astype(np.int64)               # retrigger veto, see cm_trigger_test
            dbc = bc[..., :, None] - bc[..., None, :]
            valid &= ~((dbc >= 1) & (dbc <= 3) & valid[..., None, :]).any(-1)
            good = valid[..., None] & np.isfinite(P) & np.isfinite(PL) & np.isfinite(M) & (adc > 0)
            r = adc - P
            rl = d["adc_low"].astype(np.float32) - PL
            hbg = hb * good
            A = np.where((slabs == 12)[None, :, None, None, None], A2c[cidc], A2[cidc])
            far = good & (hb == 0) & ~(np.einsum("...sc,...cd->...sd", hbg, A) > 0)
            Bm = np.where(far.sum(-1) >= 3, np.nanmedian(np.where(far, r, np.float32(np.nan)), axis=-1), np.nan)
            qadc = np.where(hbg > 0, rl, 0) / K_ANCHOR                           # HG ADC
            Qtot = qadc.sum(-1)
            B = np.where(np.isfinite(Bm), Bm, B_QUIET + B_SLOPE * Qtot)          # [a, slab, chip, sca]
            Blg = np.where(far.sum(-1) >= 3, np.nanmedian(np.where(far, rl, np.float32(np.nan)), axis=-1), np.nan)
            bmeas += [np.isfinite(Bm[Qtot > 0]).sum(), (Qtot > 0).sum()]
            qm = qadc / M
            Q = qm.sum(-1)[..., None] - qm                                       # MIP, without the channel itself
            use = good & (good.sum(-1) > 8)[..., None]
            kq = np.clip(np.searchsorted(Q_EDGES, Q[use], side="right") - 1, 0, nq - 1)
            rc = r - B[..., None]
            bit = hb[use] > 0
            for key, x, sel in (("raw_all", r / M, None), ("cor_all", rc / M, None),
                                ("cor_bit", rc / M, bit), ("cor_nobit", rc / M, ~bit)):
                xv = x[use]
                i = np.floor((xv - S_BINS[0]) / (S_BINS[1] - S_BINS[0])).astype(int)
                ok = (i >= 0) & (i < ns)
                if sel is not None:
                    ok &= sel
                h[key] += np.bincount((kq * ns + i)[ok], minlength=nq * ns)
            lin = bit & (rc[use] > 150) & (rc[use] < 800)
            blg_u = np.broadcast_to(Blg[..., None], r.shape)[use]
            lin &= np.isfinite(blg_u)   # both ratios on the chip-SCAs where the LG baseline is measured
            np.add.at(blg, kq[lin], np.stack([blg_u[lin], np.ones(lin.sum())], axis=1))
            for key, num in (("ratio", rl[use][lin] - blg_u[lin]), ("ratio_raw", rl[use][lin])):
                i = np.floor((num / rc[use][lin] - R_BINS[0]) / (R_BINS[1] - R_BINS[0])).astype(int)
                ok = (i >= 0) & (i < nr)
                hr[key] += np.bincount((kq[lin] * nr + i)[ok], minlength=nq * nr)
            n_acq += len(cid)
        print(f"[{a.run}] {os.path.basename(path)}: {n_acq} acquisitions", flush=True)
    np.savez_compressed(f"{a.outdir}/cm_mip_peak_{a.run}.npz", q_edges=Q_EDGES, s_bins=S_BINS, r_bins=R_BINS,
                        n_acq=n_acq, b_measured_frac=bmeas[0] / max(bmeas[1], 1),
                        b_lg=blg[:, 0] / np.maximum(blg[:, 1], 1),
                        **{k: v.reshape(nq, ns) for k, v in h.items()}, **{k: v.reshape(nq, nr) for k, v in hr.items()})
    print(f"[{a.run}] done: {n_acq} acquisitions, baseline measured in {bmeas[0] / max(bmeas[1], 1):.2f} "
          f"of the chip-SCAs with charge")


if __name__ == "__main__":
    main()

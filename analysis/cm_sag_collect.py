#!/usr/bin/env python
"""
Records for the baseline-sag model B = k_chip * Q_chip + f(Q_slab): per chip-SCA, the baseline of its signal-free
channels together with the charge of its own chip and of the rest of its slab and of the detector at the same time.

Same selection as cm_trigger_test.py (which only kept Q_chip): r = adc_high - pedestal (EOS Calib/MuonCalib th210,
same SCA), B = median r of the non-fired channels >= 2 pads from any fired one, charge of a fired channel = low-gain
residual / K_ANCHOR (HG ADC), SCA columns that follow another column of the same chip within 1-3 BCIDs vetoed.
Charges at "the same time" are summed over the valid columns whose BCID is within +-1 of the record's column:
  q_chip   charge of the chip-SCA itself (all fired channels)
  q_slab   charge of the other chips of the same slab
  q_det    charge of the other slabs
and, over the whole acquisition (all BCIDs), to test a slow baseline recovery after earlier events:
  q_slab_acq_prev  charge of the slab in the EARLIER columns of the acquisition (BCID < record's BCID - 1)
  q_slab_acq       charge of the slab in the whole acquisition, minus the same-time window
  n_prev_slab      number of earlier distinct BCIDs of the slab in the acquisition

    python analysis/cm_sag_collect.py <outdir> --run e286 --chunks 10-17

Writes cm_sag_<run>.npz: rec [N, 12] float32 with rec_cols
(slab, chip, sca, nbits, q_chip, q_slab, q_det, base_far, n_far, q_slab_acq_prev, q_slab_acq, n_prev_slab), and n_acq.
"""
import argparse
import os
import sys

import numpy as np
import uproot

sys.path.insert(0, __file__.rsplit("/siwecal_k4sim/", 1)[0] + "/siwecal-tb2026/diagnostics")
from cm_baseline_collect import load_ped, adjacency, run_name, CONV, CAL, MAP, MAP_COB, K_ANCHOR  # noqa: E402

NBC = 4096 + 2          # 12-bit BCID, one guard bin on each side for the +-1 window


def window_sum(hist):
    """hist[..., NBC] -> sum over bcid-1 .. bcid+1 (guard bins at 0 and NBC-1)."""
    out = hist.copy()
    out[..., 1:] += hist[..., :-1]
    out[..., :-1] += hist[..., 1:]
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("outdir")
    ap.add_argument("--run", required=True)
    ap.add_argument("--chunks", default="10-17")
    ap.add_argument("--max-acq", type=int, default=0)
    a = ap.parse_args()
    os.makedirs(a.outdir, exist_ok=True)
    ped_hg = load_ped(f"{CAL}/Pedestal_TB2026CERN_run_000th210_highgain.txt")
    ped_lg = load_ped(f"{CAL}/Pedestal_TB2026CERN_run_000th210_lowgain.txt")
    _, A2 = adjacency(MAP)
    _, A2c = adjacency(MAP_COB)
    c0, c1 = (int(x) for x in a.chunks.split("-"))
    files = [f"{CONV}/{run_name(a.run)}/chunks/chunk_{i:04d}.root" for i in range(c0, c1 + 1)]
    files = [f for f in files if os.path.exists(f)]
    slabs = np.arange(15)
    recs, n_acq = [], 0
    for path in files:
        for d in uproot.open(path)["siwecaldecoded"].iterate(
                ["chipid", "nColumns", "bcid", "adc_high", "adc_low", "hitbit_high"], step_size=20, library="np"):
            if a.max_acq and n_acq >= a.max_acq:
                break
            cid = d["chipid"]
            cidc = np.clip(cid, 0, 15)
            L = slabs[None, :, None]
            P, PL = ped_hg[L, cidc], ped_lg[L, cidc]                   # [a, slab, chip, sca, ch]
            adc = d["adc_high"].astype(np.float32)
            hb = d["hitbit_high"].astype(np.float32)
            valid = (np.arange(15)[None, None, None, :] < d["nColumns"][..., None]) & (cid >= 0)[..., None]
            bc = np.where(valid, d["bcid"], -999).astype(np.int64)
            dbc = bc[..., :, None] - bc[..., None, :]
            prev = ((dbc >= 1) & (dbc <= 3) & valid[..., None, :]).any(-1)
            valid &= ~prev
            good = valid[..., None] & np.isfinite(P) & np.isfinite(PL) & (adc > 0)
            r = adc - P
            rl = d["adc_low"].astype(np.float32) - PL
            hbg = hb * good
            A = np.where((slabs == 12)[None, :, None, None, None], A2c[cidc], A2[cidc])
            near2 = np.einsum("...sc,...cd->...sd", hbg, A) > 0
            far = good & (hb == 0) & ~near2
            B = np.nanmedian(np.where(far, r, np.float32(np.nan)), axis=-1)            # [a, slab, chip, sca]
            qcol = (np.where(hbg > 0, rl, 0) / K_ANCHOR).sum(-1)                         # [a, slab, chip, sca]
            # charge per (acquisition, slab, bcid), then the +-1 window
            na = len(cid)
            ia, isl, ich, isc = np.nonzero(valid)
            b1 = d["bcid"][ia, isl, ich, isc].astype(np.int64) + 1
            H = np.zeros((na, 15, NBC), np.float64)
            np.add.at(H, (ia, isl, b1), qcol[ia, isl, ich, isc])
            W = window_sum(H)                                                            # [a, slab, bcid]
            Wdet = W.sum(1)                                                              # [a, bcid]
            Cum = np.cumsum(H, axis=-1)                                                  # charge up to each bcid
            Tot = Cum[..., -1]                                                           # [a, slab]
            Occ = np.cumsum(H > 0, axis=-1)                                              # distinct bcids so far
            sel = valid & np.isfinite(B) & (good.sum(-1) > 8)
            ia, isl, ich, isc = np.nonzero(sel)
            b1 = d["bcid"][ia, isl, ich, isc].astype(np.int64) + 1
            q_chip = qcol[ia, isl, ich, isc]
            q_slab_all = W[ia, isl, b1]
            q_slab = q_slab_all - q_chip                       # other chips of the slab (own column removed)
            q_det = Wdet[ia, b1] - q_slab_all                  # other slabs
            q_prev = Cum[ia, isl, np.maximum(b1 - 2, 0)]       # strictly before the +-1 window
            q_acq = Tot[ia, isl] - q_slab_all
            n_prev = Occ[ia, isl, np.maximum(b1 - 2, 0)]
            recs.append(np.stack([slabs[isl], cid[ia, isl, ich], isc, hbg.sum(-1)[ia, isl, ich, isc], q_chip, q_slab,
                                  q_det, B[ia, isl, ich, isc], far.sum(-1)[ia, isl, ich, isc], q_prev, q_acq, n_prev],
                                 axis=1).astype(np.float32))
            n_acq += na
        print(f"[{a.run}] {os.path.basename(path)}: {n_acq} acquisitions", flush=True)
    rec = np.concatenate(recs)
    np.savez_compressed(f"{a.outdir}/cm_sag_{a.run}.npz", rec=rec, n_acq=n_acq,
                        rec_cols=np.array(["slab", "chip", "sca", "nbits", "q_chip", "q_slab", "q_det", "base_far",
                                           "n_far", "q_slab_acq_prev", "q_slab_acq", "n_prev_slab"]))
    print(f"[{a.run}] {len(rec)} chip-SCA records from {n_acq} acquisitions -> {a.outdir}/cm_sag_{a.run}.npz")


if __name__ == "__main__":
    main()

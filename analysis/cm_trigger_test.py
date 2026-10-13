#!/usr/bin/env python
"""
Does the chip baseline drop also reach the trigger? Hit-bit turn-on against the high-gain ADC, by chip-SCA charge.

In a busy chip-SCA the baseline of the channels without a hit bit sits below the pedestal (diagnostics/
cm_baseline_collect.py in siwecal-tb2026: ~-1 ADC per 1000 ADC of charge). The ADC sample of a channel is then
pedestal + baseline + signal. Two cases, told apart by where the hit-bit turn-on sits:

  trigger NOT affected (fast shaper keeps its baseline): P(bit) depends on the signal = r - baseline, so the turn-on
        is fixed in the CORRECTED residual r - baseline and moves down in the raw residual r = adc - pedestal.
  trigger affected by the same drop: a larger signal is needed to fire, so the turn-on is fixed in the RAW residual
        and moves up in the corrected one -> hits are lost at the discriminator in busy chips.

For every channel of every valid chip-SCA: r = adc_high - pedestal (EOS Calib/MuonCalib th210, same SCA), its hit
bit, and the chip-SCA charge Q (fired channels, low gain / k, without the channel itself) and baseline B (median r of
the non-fired channels >= 2 pads from any fired one). P(bit | r) is fitted per Q bin with
f0 + (p - f0) * (1 + erf((r - mu) / (sqrt2 sigma))) / 2, in r and in r - B.

    python analysis/cm_trigger_test.py <outdir> --run e291 --chunks 10-17 [--slab-chip 6:11]

Writes cm_trigger_<run>[_<slab>_<chip>].npz (histograms) and .txt (fits).
"""
import argparse
import os
import sys

import numpy as np
import uproot
from scipy.optimize import curve_fit
from scipy.special import erf

sys.path.insert(0, __file__.rsplit("/siwecal_k4sim/", 1)[0] + "/siwecal-tb2026/diagnostics")
from cm_baseline_collect import load_ped, adjacency, run_name, CONV, CAL, MAP, MAP_COB, K_ANCHOR  # noqa: E402

Q_EDGES = np.array([0, 200, 500, 1000, 2000, 4000, 8000, 1e9])
R_BINS = np.arange(-20, 60.001, 1.0)


def turnon(x, f0, p, mu, sg):
    return f0 + (p - f0) * 0.5 * (1 + erf((x - mu) / (np.sqrt(2) * sg)))


def fit(h_bit, h_all):
    ctr = 0.5 * (R_BINS[1:] + R_BINS[:-1])
    ok = h_all > 30
    if ok.sum() < 15:
        return [np.nan] * 4
    eff = h_bit[ok] / h_all[ok]
    err = np.sqrt(np.clip(eff * (1 - eff), 1e-4, None) / h_all[ok])
    try:
        p, _ = curve_fit(turnon, ctr[ok], eff, p0=[0.01, 0.95, 16, 3.5], sigma=err,
                         bounds=([0, 0.3, 0, 0.5], [0.3, 1.0, 50, 15]))
        return list(p)
    except RuntimeError:
        return [np.nan] * 4


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("outdir")
    ap.add_argument("--run", required=True)
    ap.add_argument("--chunks", default="10-17")
    ap.add_argument("--slab-chip", default=None, help="restrict to one slab:chipid")
    ap.add_argument("--no-veto", dest="veto_retrigger", action="store_false",
                    help="keep the SCA columns that follow another one within 1-3 BCIDs")
    ap.add_argument("--tag", default="", help="suffix for the output names")
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
    nq, nr = len(Q_EDGES) - 1, len(R_BINS) - 1
    NCH = 15 * 16 * 64
    h = {k: np.zeros(NCH * nq * nr, np.int32) for k in ("raw_all", "raw_bit", "cor_all", "cor_bit")}
    bsum, bn = np.zeros(nq), np.zeros(nq)
    recs = []
    n_acq = 0
    for path in files:
        for d in uproot.open(path)["siwecaldecoded"].iterate(
                ["chipid", "nColumns", "bcid", "adc_high", "adc_low", "hitbit_high"], step_size=20, library="np"):
            if a.max_acq and n_acq >= a.max_acq:
                break
            slabs = np.arange(15)
            if a.slab_chip:   # keep only that slab before the heavy part
                s0 = int(a.slab_chip.split(":")[0])
                d = {k: v[:, s0:s0 + 1] for k, v in d.items()}
                slabs = slabs[s0:s0 + 1]
            cid = d["chipid"]
            cidc = np.clip(cid, 0, 15)
            L = slabs[None, :, None]
            P, PL = ped_hg[L, cidc], ped_lg[L, cidc]                   # [a, slab, chip, sca, ch]
            adc = d["adc_high"].astype(np.float32)
            hb = d["hitbit_high"].astype(np.float32)
            valid = (np.arange(15)[None, None, None, :] < d["nColumns"][..., None]) & (cid >= 0)[..., None]
            if a.slab_chip:
                s, c = (int(x) for x in a.slab_chip.split(":"))
                valid &= ((slabs == s)[None, :, None] & (cid == c))[..., None]
            if a.veto_retrigger:
                # drop SCA columns that follow another column of the same chip within 1-3 BCIDs: there the channels
                # that fired in the earlier column still show their signal, without a hit bit
                bc = np.where(valid, d["bcid"], -999).astype(np.int64)
                dbc = bc[..., :, None] - bc[..., None, :]                                # [a, slab, chip, sca, sca']
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
            qch = np.where(hbg > 0, rl, 0) / K_ANCHOR
            Q = qch.sum(-1)[..., None] - qch                                         # without the channel itself
            use = good & np.isfinite(B)[..., None] & (good.sum(-1) > 8)[..., None]
            kq = np.clip(np.searchsorted(Q_EDGES, Q, side="right") - 1, 0, nq - 1)[use]
            rr, bb, bit = r[use], np.broadcast_to(B[..., None], r.shape)[use], hb[use] > 0
            chan = (slabs[None, :, None, None, None] * 16 + cidc[..., None, None]) * 64 + np.arange(64)
            chan = np.broadcast_to(chan, r.shape)[use]
            for tag, x in (("raw", rr), ("cor", rr - bb)):
                ir = np.floor((x - R_BINS[0]) / (R_BINS[1] - R_BINS[0])).astype(int)
                ok = (ir >= 0) & (ir < nr)
                flat = (chan * nq + kq) * nr + ir
                h[f"{tag}_all"] += np.bincount(flat[ok], minlength=NCH * nq * nr).astype(np.int32)
                h[f"{tag}_bit"] += np.bincount(flat[ok & bit], minlength=NCH * nq * nr).astype(np.int32)
            np.add.at(bsum, kq[~bit], bb[~bit]); np.add.at(bn, kq[~bit], 1)
            # per chip-SCA record: slab, chip id, hit bits, charge, baseline, number of far channels
            sel = np.isfinite(B) & (good.sum(-1) > 8)
            ia, isl, ich, isc = np.nonzero(sel)
            recs.append(np.stack([slabs[isl], cid[ia, isl, ich], hbg.sum(-1)[sel], qch.sum(-1)[sel], B[sel],
                                  far.sum(-1)[sel]], axis=1).astype(np.float32))
            n_acq += len(cid)
        print(f"[{a.run}] {os.path.basename(path)}: {n_acq} acquisitions", flush=True)
    tag = a.run + (f"_{a.slab_chip.replace(':', '_')}" if a.slab_chip else "") + a.tag
    h = {k: v.reshape(NCH, nq, nr) for k, v in h.items()}
    # channels that do not trigger on a clear signal (masked or dead discriminator): bit in < 50 % of r >= 35 ADC
    hi = R_BINS[:-1] >= 35
    n_hi = h["raw_all"][:, :, hi].sum((1, 2)); b_hi = h["raw_bit"][:, :, hi].sum((1, 2))
    masked = (n_hi >= 20) & (b_hi < 0.5 * n_hi)
    np.savez_compressed(f"{a.outdir}/cm_trigger_{tag}.npz", q_edges=Q_EDGES, r_bins=R_BINS, base=bsum / np.maximum(bn, 1),
                        n_acq=n_acq, masked=masked, rec=np.concatenate(recs),
                        rec_cols=np.array(["slab", "chip", "nbits", "q", "base_far", "n_far"]), **{k: v[~masked].sum(0) for k, v in h.items()},
                        **{f"all_{k}": v.sum(0) for k, v in h.items()})
    h = {k: v[~masked].sum(0) for k, v in h.items()}
    lines = [f"== {tag}: {n_acq} acquisitions. Hit-bit turn-on per chip-SCA charge bin (HG ADC, without the channel); "
             f"{masked.sum()} non-triggering channels excluded (of {(n_hi >= 20).sum()} with >= 20 samples at r >= 35)",
             f"{'Q bin':>13} {'<B>':>6} | {'raw mu':>7} {'sigma':>6} {'plat':>5} | {'cor mu':>7} {'sigma':>6} {'plat':>5}"]
    for k in range(nq):
        fr = fit(h["raw_bit"][k], h["raw_all"][k])
        fc = fit(h["cor_bit"][k], h["cor_all"][k])
        lines.append(f"{Q_EDGES[k]:6.0f}-{Q_EDGES[k + 1]:6.0f} {bsum[k] / max(bn[k], 1):+6.2f} | "
                     f"{fr[2]:7.2f} {fr[3]:6.2f} {fr[1]:5.2f} | {fc[2]:7.2f} {fc[3]:6.2f} {fc[1]:5.2f}")
    txt = "\n".join(lines)
    print(txt)
    open(f"{a.outdir}/cm_trigger_{tag}.txt", "w").write(txt + "\n")


if __name__ == "__main__":
    main()

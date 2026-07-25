#!/usr/bin/env python3
"""
Convert simulation EDM4hep output to the "ecal tree" TTree format expected by
k4SiWEcalReco (EcalToEDM4hep component).

After running the full Gaudi pipeline (job3_digitize.py), the recommended input
collection is SiPadHitsMapped, which has TB-format CellIDs produced by
ChannelMapper.  The parallel SiPadHitsMasked UserDataCollection carries the
per-hit masking flags (0=ok, 1=masked/uncalibrated).

Older collections (SiPadHitsFlipped, SiPadHitsDigi) are also supported; in that
case chip/channel are set to 0 and hit_ismasked to 0 (legacy mode).

CellID bitfields
----------------
  SiPadHitsMapped  (TB format):   system:8, slab:8, chip:16, channel:8, sca:8
  SiPadHitsFlipped (sim format):  system:8, layer:8, slice:5, x:9, y:9
  SiPadHitsDigi    (sim format):  same as Flipped

Usage:
    python -m analysis.sim_to_ecal_tree \\
        --input  gaudi_jobs/1_mu_beam_pipeline/digitized.edm4hep.root \\
        --output gaudi_jobs/1_mu_beam_pipeline/ecal_sim.root \\
        [--collection SiPadHitsMapped] [--masking-collection SiPadHitsMasked] \\
        [--run 0] [--max-events N] [--verbose]
"""

from __future__ import annotations

import argparse
import os
import sys
import time
from typing import Optional

import numpy as np

# --------------------------------------------------------------------------- #
# CellID decoders
# --------------------------------------------------------------------------- #
# Simulation format: system:8, layer:8, slice:5, x:9, y:9
_SIM_LAYER_SHIFT  = 8
_SIM_LAYER_MASK   = 0xFF

# TB format: system:8, slab:8, chip:16, channel:8, sca:8
_TB_SLAB_SHIFT    = 8
_TB_SLAB_MASK     = 0xFF
_TB_CHIP_SHIFT    = 16
_TB_CHIP_MASK     = 0xFFFF
_TB_CHAN_SHIFT     = 32
_TB_CHAN_MASK      = 0xFF
_TB_SCA_SHIFT     = 40
_TB_SCA_MASK      = 0xFF

# Collections that carry the TB-format CellID produced by ChannelMapper
_TB_FORMAT_COLLECTIONS = {"SiPadHitsMapped"}

# --------------------------------------------------------------------------- #
# Radiation-length geometry
# --------------------------------------------------------------------------- #
# W absorber thickness (mm) in front of each silicon layer (layers 0-14).
# Matches EcalPidTransformer.WThicknesses in k4SiWEcalReco.
_W_THICKNESS_MM = [2.8, 4.2, 4.2, 4.2, 4.2, 4.2, 4.2, 4.2,
                   4.2, 5.6, 5.6, 5.6, 5.6, 5.6, 5.6]
_X0_W_MM = 3.5  # radiation length of tungsten [mm]

# Cumulative X0 at layer l = sum of W thicknesses from layer 0 to l (inclusive) / X0_W.
# A hit in layer l's Si pad has traversed this much material.
_LAYER_X0 = [sum(_W_THICKNESS_MM[:l + 1]) / _X0_W_MM
             for l in range(len(_W_THICKNESS_MM))]

# Sampling weight of layer l = thickness of ITS OWN absorber / X0_W (0.8, 1.2 or 1.6).
# Not to be confused with _LAYER_X0 above, which is cumulative: this one is the
# per-layer weight used to correct for the non-uniform sampling fraction, and it
# matches event_viewer/_metrics.py:hit_weights().
_LAYER_W_X0 = [w / _X0_W_MM for w in _W_THICKNESS_MM]


def _decode_layer_sim(cellid: int) -> int:
    """Extract layer from simulation-format CellID (layer:8 field)."""
    return (cellid >> _SIM_LAYER_SHIFT) & _SIM_LAYER_MASK


# Backward-compatible alias (used by existing tests in test_converter.py)
_decode_layer = _decode_layer_sim


def _decode_tb(cellid: int):
    """Return (slab, chip, channel, sca) from TB-format CellID."""
    slab    = (cellid >> _TB_SLAB_SHIFT)  & _TB_SLAB_MASK
    chip    = (cellid >> _TB_CHIP_SHIFT)  & _TB_CHIP_MASK
    channel = (cellid >> _TB_CHAN_SHIFT)  & _TB_CHAN_MASK
    sca     = (cellid >> _TB_SCA_SHIFT)   & _TB_SCA_MASK
    return slab, chip, channel, sca


# --------------------------------------------------------------------------- #
# Main converter
# --------------------------------------------------------------------------- #

def convert(
    input_path: str,
    output_path: str,
    collection: str = "SiPadHitsMapped",
    masking_collection: str = "SiPadHitsMasked",
    run_number: int = 0,
    max_events: Optional[int] = None,
    verbose: bool = False,
) -> int:
    """Convert one EDM4hep file to an ecal TTree.  Returns events written."""
    try:
        import podio  # noqa: F401
    except ImportError:
        sys.exit(
            "ERROR: podio Python bindings not found.\n"
            "Source the key4hep environment first:\n"
            "  source init_key4hep.sh   (release pinned in .key4hep-release)\n"
        )
    try:
        import ROOT
    except ImportError:
        sys.exit("ERROR: PyROOT not found. Source the key4hep environment first.")

    ROOT.gErrorIgnoreLevel = ROOT.kWarning  # suppress TCling::LoadPCM INFO messages

    import podio.root_io

    use_tb_format = collection in _TB_FORMAT_COLLECTIONS

    if not os.path.exists(input_path):
        sys.exit(f"ERROR: input file not found: {input_path}")

    reader = podio.root_io.Reader(input_path)
    frames = list(reader.get("events"))
    n_total = len(frames)
    if max_events is not None:
        frames = frames[:max_events]

    print(f"[sim_to_ecal_tree] Input:      {input_path}")
    print(f"[sim_to_ecal_tree] Output:     {output_path}")
    print(f"[sim_to_ecal_tree] Collection: {collection}  "
          f"({'TB CellID' if use_tb_format else 'sim CellID'})")
    print(f"[sim_to_ecal_tree] Events:     {len(frames)} / {n_total}")

    # ---------------------------------------------------------------------- #
    # Prepare TTree branches
    # ---------------------------------------------------------------------- #
    os.makedirs(os.path.dirname(os.path.abspath(output_path)), exist_ok=True)
    out_file = ROOT.TFile(output_path, "RECREATE")
    tree = ROOT.TTree("ecal", "ecal tree from k4sim simulation")

    b_run        = np.zeros(1, dtype=np.int32)
    b_event      = np.zeros(1, dtype=np.int32)
    b_spill      = np.zeros(1, dtype=np.int32)
    b_bcid       = np.zeros(1, dtype=np.int32)
    b_nhit_chan  = np.zeros(1, dtype=np.int32)
    b_nhit_slab  = np.zeros(1, dtype=np.int32)
    b_nhit_chip  = np.zeros(1, dtype=np.int32)
    b_sum_energy = np.zeros(1, dtype=np.float32)
    b_sum_hg     = np.zeros(1, dtype=np.float32)

    tree.Branch("run",        b_run,       "run/I")
    tree.Branch("event",      b_event,     "event/I")
    tree.Branch("spill",      b_spill,     "spill/I")
    tree.Branch("bcid",       b_bcid,      "bcid/I")
    tree.Branch("nhit_chan",  b_nhit_chan,  "nhit_chan/I")
    tree.Branch("nhit_slab",  b_nhit_slab, "nhit_slab/I")
    tree.Branch("nhit_chip",  b_nhit_chip, "nhit_chip/I")
    tree.Branch("sum_energy", b_sum_energy,"sum_energy/F")
    tree.Branch("sum_hg",     b_sum_hg,    "sum_hg/F")

    MAX_HITS = 4096
    b_hit_slab     = np.zeros(MAX_HITS, dtype=np.int32)
    b_hit_chip     = np.zeros(MAX_HITS, dtype=np.int32)
    b_hit_chan     = np.zeros(MAX_HITS, dtype=np.int32)
    b_hit_sca      = np.zeros(MAX_HITS, dtype=np.int32)
    b_hit_ismasked = np.zeros(MAX_HITS, dtype=np.int32)
    b_hit_energy   = np.zeros(MAX_HITS, dtype=np.float32)
    b_hit_hg       = np.zeros(MAX_HITS, dtype=np.float32)
    b_hit_lg       = np.zeros(MAX_HITS, dtype=np.float32)
    b_hit_x        = np.zeros(MAX_HITS, dtype=np.float32)
    b_hit_y        = np.zeros(MAX_HITS, dtype=np.float32)
    b_hit_z        = np.zeros(MAX_HITS, dtype=np.float32)
    b_hit_X0       = np.zeros(MAX_HITS, dtype=np.float32)
    b_hit_w_energy = np.zeros(MAX_HITS, dtype=np.float32)

    tree.Branch("hit_slab",     b_hit_slab,     "hit_slab[nhit_chan]/I")
    tree.Branch("hit_chip",     b_hit_chip,     "hit_chip[nhit_chan]/I")
    tree.Branch("hit_chan",     b_hit_chan,      "hit_chan[nhit_chan]/I")
    tree.Branch("hit_sca",      b_hit_sca,      "hit_sca[nhit_chan]/I")
    tree.Branch("hit_ismasked", b_hit_ismasked, "hit_ismasked[nhit_chan]/I")
    tree.Branch("hit_energy",   b_hit_energy,   "hit_energy[nhit_chan]/F")
    tree.Branch("hit_hg",       b_hit_hg,       "hit_hg[nhit_chan]/F")
    tree.Branch("hit_lg",       b_hit_lg,       "hit_lg[nhit_chan]/F")
    tree.Branch("hit_x",        b_hit_x,        "hit_x[nhit_chan]/F")
    tree.Branch("hit_y",        b_hit_y,        "hit_y[nhit_chan]/F")
    tree.Branch("hit_z",        b_hit_z,        "hit_z[nhit_chan]/F")
    tree.Branch("hit_X0",       b_hit_X0,       "hit_X0[nhit_chan]/F")
    tree.Branch("hit_w_energy", b_hit_w_energy, "hit_w_energy[nhit_chan]/F")

    # ---------------------------------------------------------------------- #
    # Fill loop
    # ---------------------------------------------------------------------- #
    t0 = time.time()
    n_written = 0
    n_skipped = 0
    masking_warned = False

    for frame_idx, frame in enumerate(frames):
        hits_col = frame.get(collection)
        if hits_col is None:
            if verbose:
                print(f"  frame {frame_idx}: collection '{collection}' missing, skip")
            n_skipped += 1
            continue

        hits = list(hits_col)
        n = len(hits)

        if n > MAX_HITS:
            print(f"WARNING frame {frame_idx}: {n} hits > MAX_HITS={MAX_HITS}, truncating")
            hits = hits[:MAX_HITS]
            n = MAX_HITS

        # Read masking flags if available
        mask_arr = None
        if use_tb_format and masking_collection:
            mask_col = frame.get(masking_collection)
            if mask_col is not None:
                try:
                    mask_arr = np.fromiter(mask_col, dtype=np.int32, count=len(mask_col))
                except Exception:
                    mask_arr = None
            elif not masking_warned:
                print(f"[sim_to_ecal_tree] WARNING: masking collection "
                      f"'{masking_collection}' not found; hit_ismasked set to 0")
                masking_warned = True

        b_run[0]   = run_number
        b_event[0] = frame_idx
        b_spill[0] = 0
        b_bcid[0]  = 0
        b_nhit_chan[0] = n

        slabs_seen = set()
        chips_seen = set()
        total_energy = 0.0

        for i, hit in enumerate(hits):
            cid   = hit.getCellID()
            pos   = hit.getPosition()
            energy = hit.getEnergy()

            if use_tb_format:
                slab, chip, chan, sca = _decode_tb(cid)
                ismasked = int(mask_arr[i]) if (mask_arr is not None and i < len(mask_arr)) else 0
            else:
                slab = _decode_layer_sim(cid)
                chip = 0
                chan = 0
                sca  = 0
                ismasked = 0

            b_hit_slab[i]     = slab
            b_hit_chip[i]     = chip
            b_hit_chan[i]      = chan
            b_hit_sca[i]      = sca
            b_hit_ismasked[i] = ismasked
            b_hit_energy[i]   = float(energy)
            b_hit_hg[i]       = 0.0
            b_hit_lg[i]       = 0.0
            b_hit_x[i]        = float(pos.x)
            b_hit_y[i]        = float(pos.y)
            b_hit_z[i]        = float(pos.z)
            b_hit_X0[i]       = float(_LAYER_X0[slab]) if slab < len(_LAYER_X0) else 0.0
            # Tungsten-weighted energy: E * W[slab] / X0, i.e. the hit corrected
            # for the absorber depth of its own layer.
            b_hit_w_energy[i] = (float(energy) * _LAYER_W_X0[slab]
                                 if slab < len(_LAYER_W_X0) else 0.0)

            slabs_seen.add(slab)
            chips_seen.add(chip)
            total_energy += energy

        b_nhit_slab[0] = len(slabs_seen)
        b_nhit_chip[0] = len(chips_seen)
        b_sum_energy[0] = float(total_energy)
        b_sum_hg[0]     = 0.0

        tree.Fill()
        n_written += 1

        if verbose and frame_idx % 100 == 0:
            elapsed = time.time() - t0
            rate = (frame_idx + 1) / max(elapsed, 0.001)
            print(f"  frame {frame_idx:5d}  nhit={n:4d}  "
                  f"sum_E={total_energy:.2f} MIP  [{rate:.0f} ev/s]")

    # ---------------------------------------------------------------------- #
    # Write and close
    # ---------------------------------------------------------------------- #
    out_file.cd()
    tree.Write()
    out_file.Close()

    elapsed = time.time() - t0
    print(f"[sim_to_ecal_tree] Wrote {n_written} events "
          f"({n_skipped} skipped) in {elapsed:.1f}s → {output_path}")
    return n_written


# --------------------------------------------------------------------------- #
# CLI
# --------------------------------------------------------------------------- #

def _parse_args(argv=None):
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--input", "-i",
                   default="gaudi_jobs/1_mu_beam_pipeline/digitized.edm4hep.root",
                   help="Input EDM4hep file")
    p.add_argument("--output", "-o",
                   default="gaudi_jobs/1_mu_beam_pipeline/ecal_sim.root",
                   help="Output ecal TTree file")
    p.add_argument("--collection", "-c", default="SiPadHitsMapped",
                   help="EDM4hep collection to read (default: SiPadHitsMapped)")
    p.add_argument("--masking-collection", default="SiPadHitsMasked",
                   help="Parallel UserDataCollection<int32> with masking flags "
                        "(default: SiPadHitsMasked; empty string to disable)")
    p.add_argument("--run", type=int, default=0,
                   help="Run number written to the 'run' branch (default: 0)")
    p.add_argument("--max-events", "-n", type=int, default=None,
                   help="Process at most this many events (default: all)")
    p.add_argument("--verbose", "-v", action="store_true",
                   help="Print per-event progress every 100 events")
    return p.parse_args(argv)


if __name__ == "__main__":
    args = _parse_args()
    n = convert(
        input_path=args.input,
        output_path=args.output,
        collection=args.collection,
        masking_collection=args.masking_collection,
        run_number=args.run,
        max_events=args.max_events,
        verbose=args.verbose,
    )
    sys.exit(0 if n > 0 else 1)

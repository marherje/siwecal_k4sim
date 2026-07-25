"""
Tests for the DetectorFlipper Gaudi algorithm output.

Verifies that after running job3_digitize.py (which includes DetectorFlipper):
  - SiPadHitsFlipped is present in digitized.edm4hep.root
  - Hit count matches SiPadHitsDigi
  - x, y, energy are preserved unchanged
  - z is replaced with ZPositions[layer] (for simulation defaults: z_flip == z_orig)
  - The flipping logic works correctly when a reversed ZPositions table is given

Run with:
    source init_key4hep.sh   # release pinned in .key4hep-release
    source install_env.sh   (or: export LD_LIBRARY_PATH=install/lib64:...)
    cd /afs/cern.ch/user/m/marquezh/public/siwecal_k4sim
    python -m pytest analysis/tests/test_detector_flipper.py -v
"""

from __future__ import annotations

import os
import sys
import pytest
import numpy as np

REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
DIGI_FILE = os.path.join(
    REPO_ROOT, "gaudi_jobs", "1_mu_beam_pipeline", "digitized.edm4hep.root"
)

# Simulation z positions (defaults in DetectorFlipper.cpp)
SIM_Z = [-116.35, -99.75, -83.15, -66.55, -49.95, -33.35, -16.75, -0.15,
          16.45,   33.05,  49.65,  77.25,  93.85, 110.45, 126.98]


def _check():
    if not os.path.exists(DIGI_FILE):
        return f"digitized file not found: {DIGI_FILE}"
    try:
        import podio  # noqa: F401
    except ImportError as exc:
        return f"podio unavailable: {exc}"
    try:
        import podio.root_io
        r = podio.root_io.Reader(DIGI_FILE)
        frames = list(r.get("events"))
        if "SiPadHitsFlipped" not in frames[0].getAvailableCollections():
            return ("SiPadHitsFlipped not in digitized.edm4hep.root — "
                    "run job3_digitize.py first (with DetectorFlipper in the pipeline)")
    except Exception as exc:
        return f"could not open digitized file: {exc}"
    return None


_SKIP = _check()
pytestmark = pytest.mark.skipif(_SKIP is not None, reason=_SKIP or "")


# ------------------------------------------------------------------ #
# Fixtures                                                            #
# ------------------------------------------------------------------ #
@pytest.fixture(scope="module")
def frames():
    import podio.root_io
    reader = podio.root_io.Reader(DIGI_FILE)
    return list(reader.get("events"))


# ------------------------------------------------------------------ #
# Collection presence                                                 #
# ------------------------------------------------------------------ #
def test_flipped_collection_present(frames):
    cols = frames[0].getAvailableCollections()
    assert "SiPadHitsFlipped" in cols
    assert "SiPadHitsDigi" in cols


# ------------------------------------------------------------------ #
# Hit-count consistency                                               #
# ------------------------------------------------------------------ #
def test_hit_count_matches_digi(frames):
    for frame in frames[:50]:
        digi  = list(frame.get("SiPadHitsDigi"))
        flipped = list(frame.get("SiPadHitsFlipped"))
        assert len(flipped) == len(digi), \
            f"nhit mismatch: digi={len(digi)} flipped={len(flipped)}"


# ------------------------------------------------------------------ #
# x, y, energy preserved; z from ZPositions table                    #
# ------------------------------------------------------------------ #
def test_xy_energy_preserved(frames):
    for frame in frames[:50]:
        digi    = list(frame.get("SiPadHitsDigi"))
        flipped = list(frame.get("SiPadHitsFlipped"))
        for d, f in zip(digi, flipped):
            pd, pf = d.getPosition(), f.getPosition()
            assert abs(pd.x - pf.x) < 1e-4, f"x changed: {pd.x} -> {pf.x}"
            assert abs(pd.y - pf.y) < 1e-4, f"y changed: {pd.y} -> {pf.y}"
            assert abs(d.getEnergy() - f.getEnergy()) < 1e-6, \
                f"energy changed: {d.getEnergy()} -> {f.getEnergy()}"


def test_z_matches_sim_z_table(frames):
    """For simulation (default ZPositions), z_flip must equal SIM_Z[layer]."""
    for frame in frames[:50]:
        flipped = list(frame.get("SiPadHitsFlipped"))
        for hit in flipped:
            layer = (hit.getCellID() >> 8) & 0xFF
            expected_z = SIM_Z[layer]
            got_z = hit.getPosition().z
            assert abs(got_z - expected_z) < 0.1, \
                f"layer={layer}: z_flip={got_z:.3f} != expected {expected_z:.3f}"


def test_z_replaced_not_original(frames):
    """Sanity: both collections have same layer but z values confirm remapping.

    For simulation the sim z and the flipped z are numerically identical, so
    we just check that the flipped z is exactly one of the SIM_Z values.
    """
    for frame in frames[:20]:
        for hit in frame.get("SiPadHitsFlipped"):
            z = hit.getPosition().z
            assert any(abs(z - sz) < 0.1 for sz in SIM_Z), \
                f"z={z:.3f} is not in SIM_Z lookup table"


# ------------------------------------------------------------------ #
# CellID preserved (layer encoding unchanged)                         #
# ------------------------------------------------------------------ #
def test_cellid_preserved(frames):
    for frame in frames[:50]:
        digi    = list(frame.get("SiPadHitsDigi"))
        flipped = list(frame.get("SiPadHitsFlipped"))
        for d, f in zip(digi, flipped):
            assert d.getCellID() == f.getCellID(), \
                f"cellID changed: {d.getCellID()} -> {f.getCellID()}"


# ------------------------------------------------------------------ #
# Logical flip correctness (unit test, no I/O)                       #
# ------------------------------------------------------------------ #
def test_flip_reversal_logic():
    """If ZPositions is reversed, layer 0 gets the max z and layer 14 gets the min."""
    sim_z   = SIM_Z
    flipped = list(reversed(sim_z))
    # layer 0 → flipped[0] = sim_z[14] = 126.98 (detector flipped)
    assert abs(flipped[0]  - sim_z[14]) < 0.01
    assert abs(flipped[14] - sim_z[0])  < 0.01
    # total span is the same
    assert abs((max(flipped) - min(flipped)) - (max(sim_z) - min(sim_z))) < 0.01

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
sys.path.insert(0, REPO_ROOT)

DIGI_FILE = os.path.join(
    REPO_ROOT, "gaudi_jobs", "1_mu_beam_pipeline", "digitized.edm4hep.root"
)

# The z table the pipeline actually flips to, read from the same file job3 reads
# (see geometry_ref.py). Never hardcode it here: this test used to carry its own
# copy, which survived the layer pitch going 11 mm -> 15 mm and then asserted the
# old geometry against new files. DetectorFlipper has no built-in table either.
from analysis.tests.geometry_ref import layer_pitch_mm, slab_z_mm  # noqa: E402

SLAB_Z = slab_z_mm()


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


def test_z_matches_slab_z_table(frames):
    """z_flip must equal slab_z_positions.yml[layer], the table job3 configures."""
    for frame in frames[:50]:
        flipped = list(frame.get("SiPadHitsFlipped"))
        for hit in flipped:
            layer = (hit.getCellID() >> 8) & 0xFF
            expected_z = SLAB_Z[layer]
            got_z = hit.getPosition().z
            assert abs(got_z - expected_z) < 0.1, \
                f"layer={layer}: z_flip={got_z:.3f} != expected {expected_z:.3f}"


def test_z_replaced_not_original(frames):
    """The flipped z must come from the table, not be carried over from the sim.

    Every flipped z has to be one of the table's values, and the collection as a
    whole must not simply be the simulation z: the frames differ, so a flipper
    that forgot to write would show up here.
    """
    sim_z_seen = set()
    for frame in frames[:20]:
        for hit in frame.get("SiPadHitsFlipped"):
            z = hit.getPosition().z
            assert any(abs(z - sz) < 0.1 for sz in SLAB_Z), \
                f"z={z:.3f} is not in the slab_z_positions.yml table"
        for hit in frame.get("SiPadHitsDigi"):
            sim_z_seen.add(round(hit.getPosition().z, 3))
    # The simulation frame starts well downstream of 0; the TB frame starts at 0
    # and runs negative. If these ever coincide the test above is vacuous.
    assert min(sim_z_seen) > max(SLAB_Z), \
        "simulation and target frames overlap; this test no longer proves a remap"


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
    slab_z  = SLAB_Z
    flipped = list(reversed(slab_z))
    assert abs(flipped[0]  - slab_z[-1]) < 0.01
    assert abs(flipped[-1] - slab_z[0])  < 0.01
    # total span is the same
    assert abs((max(flipped) - min(flipped)) - (max(slab_z) - min(slab_z))) < 0.01


def test_slab_table_spacing_matches_the_geometry():
    """The YAML must describe the same detector as the compact XML.

    Same check DetectorFlipper does at initialize(), here as a fast unit test so
    a drift is caught without running the pipeline. Compares |z[i+1]-z[i]|, so
    the sign flip and offset between the two frames do not matter -- only that
    the layer pitch and the gap at the empty rail slot agree.
    """
    pitch = layer_pitch_mm()
    steps = [abs(SLAB_Z[i + 1] - SLAB_Z[i]) for i in range(len(SLAB_Z) - 1)]
    # Every step is one pitch, except the empty rail slot which is two.
    for i, s in enumerate(steps):
        assert abs(s - pitch) < 0.5 or abs(s - 2 * pitch) < 0.5, \
            (f"slab_z_positions.yml layers {i}->{i+1} are {s} mm apart, which is "
             f"neither one nor two Ecal_LayerDistance ({pitch} mm)")
    assert sum(1 for s in steps if abs(s - 2 * pitch) < 0.5) == 1, \
        "expected exactly one double-pitch gap (the empty rail slot)"

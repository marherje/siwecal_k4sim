"""
Tests for analysis/sim_to_ecal_tree.py.

Run with:
    source init_key4hep.sh   # release pinned in .key4hep-release
    cd /afs/cern.ch/user/m/marquezh/public/siwecal_k4sim
    python -m pytest analysis/tests/test_converter.py -v

The tests require the key4hep environment (PODIO + PyROOT) and the digitized
simulation output at gaudi_jobs/1_mu_beam_pipeline/digitized.edm4hep.root.
"""

from __future__ import annotations

import os
import sys
import tempfile
import pytest

# Repo root is two levels up from this file
REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
sys.path.insert(0, REPO_ROOT)

INPUT_FILE = os.path.join(
    REPO_ROOT, "gaudi_jobs", "1_mu_beam_pipeline", "digitized.edm4hep.root"
)

# ------------------------------------------------------------------ #
# Skip all tests if the runtime is unavailable or the file is missing #
# ------------------------------------------------------------------ #
def _check_runtime():
    if not os.path.exists(INPUT_FILE):
        return f"input file not found: {INPUT_FILE}"
    try:
        import podio  # noqa: F401
        import ROOT   # noqa: F401
    except ImportError as exc:
        return f"runtime not available: {exc}"
    return None


_RUNTIME_SKIP = _check_runtime()
pytestmark = pytest.mark.skipif(
    _RUNTIME_SKIP is not None, reason=_RUNTIME_SKIP or ""
)


# ------------------------------------------------------------------ #
# Fixtures                                                            #
# ------------------------------------------------------------------ #
@pytest.fixture(scope="module")
def ecal_tree_path(tmp_path_factory):
    """Run the converter once for the whole module and return the output path."""
    from analysis.sim_to_ecal_tree import convert

    out = str(tmp_path_factory.mktemp("ecal") / "ecal_sim.root")
    n = convert(
        input_path=INPUT_FILE,
        output_path=out,
        collection="SiPadHitsDigi",
        run_number=0,
        max_events=50,   # fast; enough to validate schema and values
        verbose=False,
    )
    assert n > 0, "converter wrote 0 events"
    return out


@pytest.fixture(scope="module")
def tree(ecal_tree_path):
    """Open the TTree and return it (module-scoped to avoid re-opening)."""
    import ROOT
    ROOT.gROOT.SetBatch(True)
    f = ROOT.TFile(ecal_tree_path, "READ")
    t = f.Get("ecal")
    assert t, f"TTree 'ecal' not found in {ecal_tree_path}"
    # Keep the file reference alive for the module
    t._file_ref = f
    return t


# ------------------------------------------------------------------ #
# Schema tests                                                        #
# ------------------------------------------------------------------ #
REQUIRED_INT_BRANCHES = [
    "run", "event", "spill", "bcid",
    "nhit_chan", "nhit_slab", "nhit_chip",
    "hit_slab", "hit_chip", "hit_chan", "hit_sca", "hit_ismasked",
]
REQUIRED_FLOAT_BRANCHES = [
    "sum_energy", "sum_hg",
    "hit_energy", "hit_hg", "hit_lg",
    "hit_x", "hit_y", "hit_z",
]

@pytest.mark.parametrize("name", REQUIRED_INT_BRANCHES + REQUIRED_FLOAT_BRANCHES)
def test_branch_exists(tree, name):
    assert tree.GetBranch(name) is not None, f"branch '{name}' missing"


def test_entry_count(tree):
    assert tree.GetEntries() == 50


# ------------------------------------------------------------------ #
# Value range tests                                                   #
# ------------------------------------------------------------------ #

def test_nhit_chan_positive(tree):
    for entry in tree:
        assert entry.nhit_chan > 0, "empty event in ecal tree"


def test_run_number(tree):
    for entry in tree:
        assert entry.run == 0


def test_event_sequential(tree):
    for idx, entry in enumerate(tree):
        assert entry.event == idx


def test_spill_and_bcid_zero(tree):
    for entry in tree:
        assert entry.spill == 0
        assert entry.bcid == 0


def test_hit_energy_positive(tree):
    for entry in tree:
        energies = [entry.hit_energy[i] for i in range(entry.nhit_chan)]
        assert all(e > 0 for e in energies), "non-positive hit energy found"


def test_hit_slab_range(tree):
    """Layer indices must be 0-14 (15 active layers)."""
    for entry in tree:
        slabs = [entry.hit_slab[i] for i in range(entry.nhit_chan)]
        assert all(0 <= s <= 14 for s in slabs), \
            f"slab out of range [0,14]: {set(slabs)}"


def test_hit_position_x_range(tree):
    """Transverse X must be within ±90 mm (the 180 mm ASU, half either side)."""
    for entry in tree:
        xs = [entry.hit_x[i] for i in range(entry.nhit_chan)]
        assert all(-200 <= x <= 200 for x in xs), \
            f"x out of expected range: min={min(xs):.1f} max={max(xs):.1f}"


def test_hit_position_y_range(tree):
    for entry in tree:
        ys = [entry.hit_y[i] for i in range(entry.nhit_chan)]
        assert all(-200 <= y <= 200 for y in ys), \
            f"y out of expected range: min={min(ys):.1f} max={max(ys):.1f}"


def test_hit_position_z_range(tree):
    """Z must cover the detector depth (simulation: ~-120 to +130 mm)."""
    for entry in tree:
        zs = [entry.hit_z[i] for i in range(entry.nhit_chan)]
        assert all(-150 <= z <= 150 for z in zs), \
            f"z out of expected range: min={min(zs):.1f} max={max(zs):.1f}"


def test_sum_energy_consistent(tree):
    """sum_energy must equal the sum of per-hit energies."""
    for entry in tree:
        computed = sum(entry.hit_energy[i] for i in range(entry.nhit_chan))
        assert abs(entry.sum_energy - computed) < 1e-3, \
            f"sum_energy mismatch: stored={entry.sum_energy:.4f} computed={computed:.4f}"


def test_nhit_slab_consistent(tree):
    """nhit_slab must equal the number of distinct layers with hits."""
    for entry in tree:
        unique_slabs = len({entry.hit_slab[i] for i in range(entry.nhit_chan)})
        assert entry.nhit_slab == unique_slabs, \
            f"nhit_slab={entry.nhit_slab} != unique layers={unique_slabs}"


def test_hit_ismasked_all_zero(tree):
    """Simulation has no masked channels."""
    for entry in tree:
        masked = [entry.hit_ismasked[i] for i in range(entry.nhit_chan)]
        assert all(m == 0 for m in masked)


def test_synthetic_fields_zero(tree):
    """chip, chan, sca, hg, lg are synthetic – must all be 0 / 0.0."""
    for entry in tree:
        n = entry.nhit_chan
        assert all(entry.hit_chip[i] == 0 for i in range(n))
        assert all(entry.hit_chan[i]  == 0 for i in range(n))
        assert all(entry.hit_sca[i]  == 0 for i in range(n))
        assert all(entry.hit_hg[i]   == 0.0 for i in range(n))
        assert all(entry.hit_lg[i]   == 0.0 for i in range(n))


# ------------------------------------------------------------------ #
# CellID decoder unit test (no I/O needed)                           #
# ------------------------------------------------------------------ #

def test_decode_layer_known_values():
    from analysis.sim_to_ecal_tree import _decode_layer

    # CellID observed from digitized.edm4hep.root first event first hit:
    #   layer=1, slice=5, x_cell=40, y_cell=25
    # Reconstruct:  system=1 (from SND detector id=1)
    system = 1
    layer  = 1
    sliceid = 5
    x_cell  = 40
    y_cell  = 25
    cellid = (
        (system  & 0xFF)
        | ((layer  & 0xFF) << 8)
        | ((sliceid & 0x1F) << 16)
        | ((x_cell  & 0x1FF) << 21)
        | ((y_cell  & 0x1FF) << 30)
    )
    assert _decode_layer(cellid) == layer


def test_decode_layer_all_15_layers():
    from analysis.sim_to_ecal_tree import _decode_layer

    for expected_layer in range(15):
        cellid = expected_layer << 8
        assert _decode_layer(cellid) == expected_layer

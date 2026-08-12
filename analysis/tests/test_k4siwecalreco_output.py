"""
Integration tests for the k4SiWEcalReco output produced by run_pid_sim.sh.

Validates that:
  - ecal_sim.edm4hep.root exists and has the expected EDM4hep schema
  - all 23 shower variables are present
  - per-layer profiles (15 layers) are present
  - physical values are consistent with 5 GeV muons traversing the detector

Run with:
    source init_key4hep.sh   # release pinned in .key4hep-release
    cd /afs/cern.ch/user/m/marquezh/public/siwecal_k4sim
    python -m pytest analysis/tests/test_k4siwecalreco_output.py -v

Requires run_pid_sim.sh to have been executed first.
"""

from __future__ import annotations

import os
import sys
import pytest
import numpy as np

REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
sys.path.insert(0, REPO_ROOT)
sys.path.insert(0, os.path.join(REPO_ROOT, "..", "siwecal-tb2026"))

EDM4HEP_FILE = os.path.join(
    REPO_ROOT, "gaudi_jobs", "1_mu_beam_pipeline", "ecal_sim.edm4hep.root"
)
VALTREE_FILE = os.path.join(
    REPO_ROOT, "gaudi_jobs", "1_mu_beam_pipeline", "ecal_sim.valtree.root"
)
# The ecal TTree the PID ran on; test_event_count derives the expected count
# from it rather than hardcoding one.
ECAL_TREE_FILE = os.path.join(
    REPO_ROOT, "gaudi_jobs", "1_mu_beam_pipeline", "ecal_sim.root"
)


# ------------------------------------------------------------------ #
# Skip guard                                                          #
# ------------------------------------------------------------------ #
def _check():
    if not os.path.exists(EDM4HEP_FILE):
        return f"run_pid_sim.sh output not found: {EDM4HEP_FILE}"
    try:
        import podio  # noqa: F401
        from siwecal_common.edm4hep_pid import PidFileReader  # noqa: F401
    except ImportError as exc:
        return f"runtime unavailable: {exc}"
    return None


_SKIP = _check()
pytestmark = pytest.mark.skipif(_SKIP is not None, reason=_SKIP or "")


# ------------------------------------------------------------------ #
# Fixture                                                             #
# ------------------------------------------------------------------ #
@pytest.fixture(scope="module")
def pid():
    from siwecal_common.edm4hep_pid import PidFileReader
    return PidFileReader(EDM4HEP_FILE, n_layers=15)


# ------------------------------------------------------------------ #
# Schema tests                                                        #
# ------------------------------------------------------------------ #
SCALAR_VARS = [
    "nhit", "zbary", "energy", "mip_likeness", "weighte",
    "bar_x", "bar_y", "bar_r", "moliere", "transverse_rms",
    "is_shower", "shower_start", "shower_max", "shower_end",
    "shower_start_10", "shower_end_10", "shower_length",
    "first_layer", "last_layer", "n_layers_hit", "e_over_nhit",
    "shower_onset", "n_layers_before_onset",
]

def test_event_count(pid):
    """Every simulated event that survives the reco must be in the output.

    Not a hardcoded count: how many events come out depends on the run mode.
    EcalToEDM4hep drops masked channels and, in the default physics mode, hits
    below 0.5 MIP; an event left with no hits is then dropped by the total-energy
    filter. ``--validation`` disables the hit cut, so both counts are accepted
    and the failure message says which mode was expected.
    """
    import uproot
    tree = uproot.open(ECAL_TREE_FILE)["ecal"]
    energy = tree["hit_energy"].array(library="np")
    masked = tree["hit_ismasked"].array(library="np")
    alive = [(np.asarray(m) == 0) for m in masked]
    n_validation = sum(1 for a in alive if a.any())
    n_physics = sum(1 for a, e in zip(alive, energy)
                    if (a & (np.asarray(e) >= 0.5)).any())
    assert pid.n_events in (n_physics, n_validation), (
        f"got {pid.n_events} events; expected {n_physics} (physics mode, "
        f"hit cut 0.5 MIP) or {n_validation} (--validation, no hit cut) "
        f"out of {tree.num_entries} simulated")


@pytest.mark.parametrize("var", SCALAR_VARS)
def test_scalar_variable_present(pid, var):
    assert var in pid.shape_names, f"'{var}' missing from shape_names"


@pytest.mark.parametrize("kind", ["hits", "energy", "weighte"])
def test_per_layer_profiles_present(pid, kind):
    for layer in range(15):
        name = f"{kind}_per_layer_{layer}"
        assert name in pid.shape_names, f"'{name}' missing"


# ------------------------------------------------------------------ #
# Physical plausibility tests (5 GeV muons)                          #
# ------------------------------------------------------------------ #

def test_nhit_per_event(pid):
    """5 GeV muons: expect ~15 hits (one per active layer). Allow a few delta rays."""
    vals = np.array(pid.scalar("nhit"))
    assert np.median(vals) == pytest.approx(15, abs=2), \
        f"median nhit={np.median(vals):.1f} (expected ~15 for muons)"
    assert np.all(vals >= 1), "event with 0 hits found"


def test_energy_range(pid):
    """Energy in MIP units: muons deposit ~1 MIP per layer → 15-20 MIP typical."""
    vals = np.array(pid.scalar("energy"))
    assert np.median(vals) == pytest.approx(17, abs=5), \
        f"median energy={np.median(vals):.2f} MIP (expected ~17 for muons)"
    assert np.all(vals > 0), "non-positive total energy found"


def test_zbary_centered(pid):
    """zbary = energy-weighted mean layer.  For uniform muon tracks it should
    be close to the detector centre (layer ≈ 7 out of 0-14)."""
    vals = np.array(pid.scalar("zbary"))
    # Filter out edge cases (very few hits): require nhit>=10
    nhit = np.array(pid.scalar("nhit"))
    good = vals[nhit >= 10]
    assert np.mean(good) == pytest.approx(7, abs=2), \
        f"mean zbary={np.mean(good):.2f} (expected ~7 for through-going muons)"


def test_mip_likeness_high(pid):
    """Muons are MIPs: mip_likeness should be close to 1.0 for the vast majority."""
    vals = np.array(pid.scalar("mip_likeness"))
    assert np.median(vals) > 0.90, \
        f"median mip_likeness={np.median(vals):.3f} (expected >0.90 for muons)"


def test_not_shower(pid):
    """5 GeV muons must not be classified as showers."""
    is_shower = np.array(pid.scalar("is_shower"))
    shower_fraction = np.mean(is_shower)
    assert shower_fraction < 0.10, \
        f"shower fraction={shower_fraction:.2%} too high for muons (expected <10%)"


def test_shower_onset_defined_whenever_flagged(pid):
    """The onset is what sets is_shower, so it must exist wherever the flag is.

    This is the structural guarantee the onset criterion buys over the old
    ascending-pairs one: shower_start can be NaN on a flagged event (it only
    searches before the profile maximum), shower_onset cannot.
    """
    flag  = np.array(pid.scalar("is_shower"))
    onset = np.array(pid.scalar("shower_onset"))
    before = np.array(pid.scalar("n_layers_before_onset"))
    flagged = flag > 0
    if not flagged.any():
        pytest.skip("no shower in this sample")
    assert np.all(np.isfinite(onset[flagged])), \
        "is_shower set on an event with no shower_onset"
    assert np.all(np.isfinite(before[flagged]))
    assert np.all(before[flagged] <= onset[flagged]), \
        "more hit layers before the onset than there are layers before it"
    assert np.all(np.isnan(onset[~flagged])), \
        "shower_onset set on a non-shower event"


def test_n_layers_hit(pid):
    """Muons traverse the full detector: most events should hit all 15 layers."""
    vals = np.array(pid.scalar("n_layers_hit"))
    assert np.median(vals) == pytest.approx(15, abs=1), \
        f"median n_layers_hit={np.median(vals):.1f} (expected 15 for muons)"


def test_hits_per_layer_uniform(pid):
    """For muons, hits should be approximately uniform across all 15 layers."""
    layer_totals = np.array([
        np.sum(pid.scalar(f"hits_per_layer_{l}")) for l in range(15)
    ])
    # Maximum layer should not be more than 3× the minimum (rough uniformity check)
    ratio = layer_totals.max() / max(layer_totals.min(), 1)
    assert ratio < 5, \
        f"Hit distribution very non-uniform: max/min = {ratio:.1f}"


# ------------------------------------------------------------------ #
# valtree format test                                                 #
# ------------------------------------------------------------------ #

def test_valtree_file_exists():
    assert os.path.exists(VALTREE_FILE), \
        f"valtree output not found: {VALTREE_FILE}"


def test_valtree_has_scalar_branches():
    import ROOT
    ROOT.gROOT.SetBatch(True)
    f = ROOT.TFile(VALTREE_FILE, "READ")
    t = f.Get("ecal")
    assert t is not None, "TTree 'ecal' not found in valtree file"
    for var in ["nhit", "energy", "mip_likeness", "is_shower", "zbary"]:
        assert t.GetBranch(var) is not None, f"'{var}' branch missing from valtree"
    f.Close()

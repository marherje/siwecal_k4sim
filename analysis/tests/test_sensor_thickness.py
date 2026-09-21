"""
The sensor thickness per layer of the simulated geometry must be the one the
per-slab description says (mappings/slab_z_positions.yml, sensor_thickness_um).

Both are hand-maintained, and the two disagreeing is exactly the kind of thing
that shows up months later as one layer 30% off: the MIP scale of a layer is
set by its silicon thickness, and MIP values come from the geometry
(gaudi_jobs/mip_extraction_pipeline) while the thickness table drives the
per-slab treatment downstream.
"""

from __future__ import annotations

import os
import sys
import xml.etree.ElementTree as ET

import pytest

REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
GEOM_DIR = os.path.join(REPO_ROOT, "simulation", "geometry")
sys.path.insert(0, REPO_ROOT)
sys.path.insert(0, GEOM_DIR)

from analysis.slab_description import DEFAULT_YAML, load_slab_description  # noqa: E402
from parse_geometry import SiWEcalGeometry  # noqa: E402


def sensitive_thickness_per_layer(compact: str, detector: str):
    """[thickness_um per layer] of the sensitive slice of every <layer> block,
    expanded by its repeat count, with the compact's constants resolved."""
    consts = SiWEcalGeometry(os.path.join(GEOM_DIR, compact)).constants
    root = ET.parse(os.path.join(GEOM_DIR, detector)).getroot()
    out = []
    for layer in root.iter("layer"):
        repeat = int(layer.get("repeat", "1"))
        sens = [s for s in layer.findall("slice") if s.get("sensitive") == "yes"]
        assert len(sens) == 1, "one sensitive slice per layer block"
        expr = sens[0].get("thickness")
        value_mm = consts[expr] if expr in consts else float(eval(expr, {"__builtins__": {}}, dict(consts, mm=1.0)))  # noqa: S307
        out.extend([round(value_mm * 1000.0)] * repeat)
    return out


@pytest.mark.parametrize("compact,detector", [
    ("SND_compact.xml", "SiPadDetector.xml"),
    ("SND_compact_beamline.xml", "SiPadDetector_beamline.xml"),
])
def test_geometry_thickness_matches_the_slab_description(compact, detector):
    slabs = load_slab_description(DEFAULT_YAML)
    geom = sensitive_thickness_per_layer(compact, detector)
    assert len(geom) == 15
    assert geom == slabs.sensor_thickness_um, (
        f"{detector}: sensitive slice per layer {geom} != slab_z_positions.yml "
        f"{slabs.sensor_thickness_um}")


def test_layer_12_has_its_own_constant():
    """Slab 12 is the chip-on-board; its sensor thickness must be switchable on
    its own without touching the other layers."""
    consts = SiWEcalGeometry(os.path.join(GEOM_DIR, "SND_compact.xml")).constants
    assert "Ecal_WaferThickness_L12" in consts
    assert "Ecal_w_slab_gap_L12" in consts
    # the gap absorbs the thickness so the pitch stays what it is
    assert consts["Ecal_w_slab_gap_L12"] + consts["Ecal_WaferThickness_L12"] == pytest.approx(
        consts["Ecal_w_slab_gap500"] + consts["Ecal_WaferThickness500"])


def test_mip_values_follow_the_thickness():
    """The extracted MIP values (mappings/mip_values_sim.yml) scale with the
    sensor thickness: the 650 um layer collects ~1.3x the charge."""
    import yaml
    path = os.path.join(REPO_ROOT, "mappings", "mip_values_sim.yml")
    if not os.path.isfile(path):
        pytest.skip("mip_values_sim.yml not produced yet (run mip_extraction_pipeline.sh)")
    values = yaml.safe_load(open(path))["mip_gev"]
    slabs = load_slab_description(DEFAULT_YAML)
    assert len(values) == 15
    ref = [v for v, t in zip(values, slabs.sensor_thickness_um) if t == 500]
    mean500 = sum(ref) / len(ref)
    for layer, (v, t) in enumerate(zip(values, slabs.sensor_thickness_um)):
        expected = mean500 * t / 500.0
        assert abs(v / expected - 1.0) < 0.05, f"layer {layer}: {v:.3e} vs {expected:.3e}"

"""Shared geometry accessors for the analysis tests.

Tests must not carry their own copy of the detector's longitudinal layout.
Several did, and they kept asserting an 11 mm layer pitch (and, older still, a
16.6 mm one) long after the geometry moved to 15 mm — failing for a reason that
had nothing to do with the code under test.

Two sources here, and only these two anywhere:

* ``layer_pitch_mm()``  the compact XML, i.e. what the simulation was built from
* ``slab_z_mm()``       mappings/slab_z_positions.yml, the canonical test-beam
                        frame, shared with the event viewer and k4SiWEcalReco,
                        and the array job3 hands to DetectorFlipper
"""

from __future__ import annotations

import os
import xml.etree.ElementTree as ET

REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
COMPACT_FILE = os.path.join(REPO_ROOT, "simulation", "geometry", "SND_compact.xml")
SLAB_Z_FILE = os.path.join(REPO_ROOT, "mappings", "slab_z_positions.yml")


def compact_constant(name: str) -> float:
    """A <constant> from the compact XML, in mm, for the simple `N*mm` forms."""
    consts = {c.get("name"): c.get("value")
              for c in ET.parse(COMPACT_FILE).getroot().iter("constant")}
    if name not in consts:
        raise KeyError(f"'{name}' not found in {COMPACT_FILE}")
    return float(consts[name].replace("*mm", "").strip())


def layer_pitch_mm() -> float:
    """Nominal layer-to-layer distance [mm] (``Ecal_LayerDistance``)."""
    return compact_constant("Ecal_LayerDistance")


def slab_z_mm() -> list:
    """Per-slab z [mm] in the test-beam frame, from slab_z_positions.yml."""
    import yaml
    with open(SLAB_Z_FILE) as f:
        return [float(z) for z in yaml.safe_load(f)["slab_z_mm"]]

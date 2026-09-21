"""
The per-slab hardware description of the test-beam stack, from
``mappings/slab_z_positions.yml`` (the same file the test-beam repo carries).

Beyond the z and tungsten tables the file may carry a technology block:
``technologies`` (per-technology defaults: sensor thickness, pad map),
``slab_technology`` (one entry per slab), ``sensor_thickness_um`` and
``threshold_dac`` (per slab, -1 = the run's own DAC).  Every part of it is
optional; without it every slab is a FEV10 with a 500 um sensor, the default
pad map and no threshold override -- which is what the code assumed before the
block existed.

This is what turns "slab 12 is the chip-on-board" from a hard-coded ``12`` in
four places into one line of YAML: job3 asks ``pad_map_overrides()`` which slabs
need their own map, the geometry test asks ``sensor_thickness_um`` what the
sensitive slice must be, and the trigger model asks ``threshold_dac`` which slabs
ran at a discriminator of their own.
"""

from __future__ import annotations

import os
from dataclasses import dataclass
from typing import Dict, List, Optional

N_SLABS = 15
DEFAULT_TECHNOLOGY = "FEV10"
DEFAULT_THICKNESS_UM = 500
DEFAULT_PAD_MAP = "fev10_rotate_chip_channel_x_y_mapping.txt"

_HERE = os.path.dirname(os.path.abspath(__file__))
DEFAULT_YAML = os.path.normpath(os.path.join(_HERE, "..", "mappings", "slab_z_positions.yml"))


@dataclass(frozen=True)
class SlabDescription:
    technology: List[str]            # per slab
    sensor_thickness_um: List[int]   # per slab
    pad_map: List[str]               # per slab, absolute path
    threshold_dac: List[int]         # per slab, -1 = none
    slab_z_mm: List[float]
    w_thickness_mm: List[float]
    has_technology_block: bool
    source: str

    def pad_map_overrides(self) -> Dict[int, str]:
        """Slabs whose pad map differs from the default one, ``{slab: path}``.

        The default is the map of the majority technology (FEV10); the
        ChannelMapper / EcalEventBuilder take exactly this dictionary as their
        per-slab override list.
        """
        default = self.default_pad_map()
        return {s: p for s, p in enumerate(self.pad_map) if p != default}

    def default_pad_map(self) -> str:
        counts: Dict[str, int] = {}
        for p in self.pad_map:
            counts[p] = counts.get(p, 0) + 1
        return max(counts, key=lambda k: (counts[k], k == self.pad_map[0]))

    def slabs_with_threshold_override(self) -> Dict[int, int]:
        return {s: d for s, d in enumerate(self.threshold_dac) if d >= 0}


def _per_slab(values: Optional[list], name: str, default, n: int = N_SLABS) -> list:
    if values is None:
        return [default] * n
    if len(values) != n:
        raise ValueError(f"{name}: expected {n} entries, got {len(values)}")
    return list(values)


def load_slab_description(path: str = DEFAULT_YAML,
                          mappings_dir: Optional[str] = None) -> SlabDescription:
    """Read the YAML; pad maps are resolved against ``mappings_dir`` (default:
    the directory of the YAML itself)."""
    import yaml

    with open(path) as fh:
        raw = yaml.safe_load(fh) or {}
    mappings_dir = mappings_dir or os.path.dirname(os.path.abspath(path))

    slab_z = _per_slab(raw.get("slab_z_mm"), "slab_z_mm", 0.0)
    w_thick = _per_slab(raw.get("w_thickness_mm"), "w_thickness_mm", 0.0)

    technologies = raw.get("technologies") or {}
    tech = _per_slab(raw.get("slab_technology"), "slab_technology", DEFAULT_TECHNOLOGY)
    has_block = "slab_technology" in raw
    for s, t in enumerate(tech):
        if t not in technologies and t != DEFAULT_TECHNOLOGY:
            raise ValueError(f"slab {s}: technology '{t}' is not declared under 'technologies'")

    def tech_default(t: str, key: str, fallback):
        entry = technologies.get(t) or {}
        return entry.get(key, fallback)

    thickness = raw.get("sensor_thickness_um")
    if thickness is None:
        thickness = [tech_default(t, "sensor_thickness_um", DEFAULT_THICKNESS_UM) for t in tech]
    thickness = [int(v) for v in _per_slab(thickness, "sensor_thickness_um", DEFAULT_THICKNESS_UM)]

    pad_map = [tech_default(t, "pad_map", DEFAULT_PAD_MAP) for t in tech]
    pad_map = [p if os.path.isabs(p) else os.path.join(mappings_dir, p) for p in pad_map]

    dac = [int(v) for v in _per_slab(raw.get("threshold_dac"), "threshold_dac", -1)]

    return SlabDescription(technology=tech, sensor_thickness_um=thickness, pad_map=pad_map,
                           threshold_dac=dac, slab_z_mm=[float(v) for v in slab_z],
                           w_thickness_mm=[float(v) for v in w_thick],
                           has_technology_block=has_block, source=os.path.abspath(path))

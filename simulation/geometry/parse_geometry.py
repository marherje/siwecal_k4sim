"""Parse SND_compact.xml <define> constants and <readout> BitField strings.

Single source of truth for the geometry numbers the Gaudi jobs need: the pixel
pitch, the transverse envelope and — above all — the readout BitField, which
the jobs used to hardcode. That hardcoded copy went stale the moment the
segmentation gained its `wafer` field, so read it from here instead.

Usage in Gaudi job files:
    import sys
    from pathlib import Path
    sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "simulation" / "geometry"))
    from parse_geometry import SiWEcalGeometry
    geo = SiWEcalGeometry()
    bitfield = geo.bitfields["SiPadHits"]
"""

import math
import xml.etree.ElementTree as ET
from pathlib import Path

_DEFAULT_COMPACT = Path(__file__).parent / "SND_compact.xml"


class SiWEcalGeometry:
    """Resolved constants and BitField strings from SND_compact.xml."""

    def __init__(self, compact_file=None):
        self._file = Path(compact_file) if compact_file else _DEFAULT_COMPACT
        self._constants: dict[str, float] = {}
        self._bitfields: dict[str, str] = {}
        self._parse()

    # ── internal ────────────────────────────────────────────────────────────

    def _parse(self):
        root = ET.parse(self._file).getroot()

        define = root.find("define")
        if define is None:
            raise ValueError(f"No <define> block in {self._file}")
        raw = {
            c.get("name"): c.get("value", "0")
            for c in define.findall("constant")
            if c.get("name")
        }
        self._constants = self._resolve(raw)

        readouts = root.find("readouts")
        if readouts is not None:
            for ro in readouts.findall("readout"):
                name = ro.get("name")
                id_el = ro.find("id")
                if name and id_el is not None and id_el.text:
                    self._bitfields[name] = id_el.text.strip()

    def _resolve(self, raw: dict) -> dict:
        # Constants reference each other in any order, so iterate until the
        # expressions stop resolving. Units are aliases evaluating to mm.
        ctx: dict = {
            "mm": 1.0,
            "cm": 10.0,
            "m": 1000.0,
            "tesla": 1.0,
            "deg": 1.0,
            "pi": math.pi,
            "floor": math.floor,
        }
        todo = dict(raw)
        for _ in range(len(todo) + 1):
            if not todo:
                break
            progress = False
            for name, expr in list(todo.items()):
                try:
                    ctx[name] = float(eval(expr, {"__builtins__": {}}, ctx))  # noqa: S307
                    del todo[name]
                    progress = True
                except Exception:
                    pass
            if not progress:
                break
        _units = {"mm", "cm", "m", "tesla", "deg", "pi", "floor"}
        return {k: v for k, v in ctx.items() if k not in _units}

    # ── raw access ──────────────────────────────────────────────────────────

    @property
    def constants(self) -> dict:
        return self._constants

    @property
    def bitfields(self) -> dict:
        return self._bitfields

    # ── SiPad (SiW-ECAL) ────────────────────────────────────────────────────

    @property
    def ecal_cell_size_x(self) -> float:
        """Pixel pitch along x [mm] — what the readout segmentation tiles with."""
        return self._constants["Ecal_CellSizeX"]

    @property
    def ecal_cell_size_y(self) -> float:
        return self._constants["Ecal_CellSizeY"]

    @property
    def ecal_dim_x(self) -> float:
        return self._constants["Ecal_dim_x"]

    @property
    def ecal_dim_y(self) -> float:
        return self._constants["Ecal_dim_y"]

    @property
    def ecal_dim_z(self) -> float:
        return self._constants["Ecal_dim_z"]

    @property
    def ecal_n_layers(self) -> int:
        return int(self._constants["Ecal_NLayers"])

    @property
    def ecal_layer_distance(self) -> float:
        return self._constants["Ecal_LayerDistance"]

    @property
    def hough_half_size(self) -> float:
        """Half-size for the Hough histogram [mm].

        The transverse envelope plus one cell of slack: a peak can only sit on
        a pad centre, so the histogram never needs to reach beyond the plane.
        """
        return 0.5 * max(self.ecal_dim_x, self.ecal_dim_y) + self.ecal_cell_size_x

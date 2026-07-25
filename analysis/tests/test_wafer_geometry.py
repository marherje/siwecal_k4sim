"""
Tests for the tiled wafer geometry of the SiPad detector.

The point of these checks is to prove that a segmented calorimeter plane with
*physical* dead regions is self-consistent, i.e. that three things agree:

  1. the geometry   -- where the silicon actually is,
  2. the readout    -- which cell a position maps to,
  3. the real detector -- the pad map used by ChannelMapper.

The layout is two-level: a plane holds NASUs x NASUs ASUs on a fixed pitch, each
ASU holds NWafers x NWafers silicon sensors butted against each other, and each
sensor is a dead rim (Ecal_WaferMarginX, no guard ring beyond it) around a
sensitive pad array subdivided by the readout segmentation. The rim volume
carries the 'wafer' physVolID; MultiSegmentation uses that ID to pick a
per-wafer offset so the x/y fields come out as GLOBAL pad indices.

The dead regions of a layer are therefore of two kinds: the rims, which are
silicon but not sensitive, and the parts of the plane no sensor covers at all.
Testing them means asking which *volume* a point lands in, not which material.

Everything below is derived from the compact file, not hardcoded, so the same
tests work for the 18 / 36 / 54 cm configurations.

One unit trap worth remembering when porting this elsewhere: DD4hep's internal
length unit is the centimetre (dd4hep::mm == 0.1), so a constant written as
"5.52*mm" in the XML comes back from constantAsDouble() as 0.552. TGeo navigation
and DDSegmentation use the same convention, so everything here stays in cm and
only converts to mm to compare against the pad map, which is written in mm.

Run with:
    source init_key4hep.sh   # release pinned in .key4hep-release
    export LD_LIBRARY_PATH=$PWD/install/lib64:$PWD/install/lib:$LD_LIBRARY_PATH
    cd /afs/cern.ch/user/m/marquezh/public/siwecal_k4sim
    python -m pytest analysis/tests/test_wafer_geometry.py -v

These checks are static: they prove the geometry and the readout agree. They do
NOT prove Geant4 actually drops hits in the dead regions -- for that, shoot muons
with simulation/geometry/dead_zone_test.py.
"""

from __future__ import annotations

import os
import pytest

REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
GEO_DIR = os.path.join(REPO_ROOT, "simulation", "geometry")
COMPACT = os.path.join(GEO_DIR, "SND_compact.xml")
PAD_MAP = os.path.join(REPO_ROOT, "mappings", "fev10_rotate_chip_channel_x_y_mapping.txt")

READOUT = "SiPadHits"
CM_TO_MM = 10.0  # geometry works in cm; the pad map is written in mm

# Name suffix of the sensitive volume built by the plugin, i.e. the pad array
# inside a sensor. Everything else a point can land in is dead.
ACTIVE_SUFFIX = "_wafer_pads"

# The pad map is written on the grid of this very geometry -- see
# mappings/regenerate_pad_maps.py -- so the two agree to the rounding of the map
# file, which holds three decimals of a millimetre.
PAD_MAP_TOL_MM = 1e-3


def _require():
    if not os.path.exists(COMPACT):
        return f"compact file not found: {COMPACT}"
    try:
        import dd4hep  # noqa: F401
        import ROOT  # noqa: F401
    except ImportError as exc:
        return f"key4hep not set up: {exc}"
    return None


pytestmark = pytest.mark.skipif(_require() is not None, reason=_require() or "")


class Layout:
    """Wafer/pad layout read straight out of the compact file, in cm."""

    def __init__(self, detector):
        c = detector.constantAsDouble
        self.pad = c("Ecal_CellSizeX")            # pixel pitch
        self.diode = c("Ecal_PadSizeX")           # active part of a pixel
        self.pad_margin = c("Ecal_PadMarginX")    # inactive part of a pixel
        self.pads_per_wafer = int(c("Ecal_NPadsPerWaferX"))
        self.active = c("Ecal_WaferActiveX")      # sensitive area of a sensor
        self.wafer = c("Ecal_WaferSizeX")         # the whole sensor
        self.margin = c("Ecal_WaferMarginX")      # its inactive rim
        self.gap = c("Ecal_WaferGapX")            # extra clearance, if any
        self.wafers_per_asu = int(c("Ecal_NWafersX"))
        self.asus = int(c("Ecal_NASUsX"))
        self.asu_pitch = c("Ecal_ASUPitchX")
        self.dim = c("Ecal_dim_x")
        self.cols = self.asus * self.wafers_per_asu          # wafer columns
        self.npads = self.cols * self.pads_per_wafer         # pads per side

    def wafer_centre(self, col: int) -> float:
        """Centre of global wafer column `col`."""
        asu, w = divmod(col, self.wafers_per_asu)
        return ((asu - (self.asus - 1) / 2.0) * self.asu_pitch
                + (w - (self.wafers_per_asu - 1) / 2.0) * (self.wafer + self.gap))

    def local(self, pad: int) -> float:
        """Centre of a pad within its sensor's *pad array*.

        That array is the sensitive volume, so it is the frame the segmentation
        sees; the sensor's rim is outside it and plays no part.
        """
        return -self.active / 2.0 + (pad % self.pads_per_wafer + 0.5) * self.pad

    def globl(self, pad: int) -> float:
        """Centre of a global pad index over the whole plane."""
        return self.wafer_centre(pad // self.pads_per_wafer) + self.local(pad)

    def wafer_id(self, ix: int, iy: int) -> int:
        """physVolID of the wafer holding global pad (ix, iy)."""
        return (iy // self.pads_per_wafer) * self.cols + ix // self.pads_per_wafer


@pytest.fixture(scope="module")
def geo():
    import dd4hep
    import ROOT

    det = dd4hep.Detector.getInstance()
    # dd4hep::Detector is a process-wide singleton and fromXML() throws if it
    # already holds a geometry. Another test module in the same pytest run may
    # have loaded this very compact file first (see test_viewer_wafer_display),
    # so only load it when nobody has.
    if det.state() != dd4hep.Detector.READY:
        det.fromXML(COMPACT)
    layout = Layout(det)

    # z of the first silicon plane, found by scanning along the beam axis
    gm = ROOT.gGeoManager
    z, z_si = 0.0, None
    while z < 60.0:
        node = gm.FindNode(layout.globl(0), layout.globl(0), z)
        if node and node.GetVolume().GetMedium().GetMaterial().GetName().startswith("Si"):
            z_si = z
            break
        z += 0.005
    assert z_si is not None, "no silicon plane found along z"
    return det, gm, layout, z_si


def _material(gm, x, y, z):
    return gm.FindNode(x, y, z).GetVolume().GetMedium().GetMaterial().GetName()


def _is_active(gm, x, y, z):
    """True if (x, y, z) lands in the sensitive pad array of some sensor.

    Material is not the question any more: a sensor's rim is silicon too, and
    asking for silicon would call the dead cross between two sensors active.
    """
    node = gm.FindNode(x, y, z)
    return bool(node) and node.GetVolume().GetName().endswith(ACTIVE_SUFFIX)


def test_every_pad_centre_is_active_silicon(geo):
    """The active area: one silicon pad per readout cell, none missing."""
    _, gm, lay, z = geo
    misses = [
        (i, j)
        for i in range(lay.npads)
        for j in range(lay.npads)
        if not (_is_active(gm, lay.globl(i), lay.globl(j), z)
                and _material(gm, lay.globl(i), lay.globl(j), z).startswith("Si"))
    ]
    assert not misses, f"{len(misses)} pad centres are not active silicon, e.g. {misses[:5]}"


def test_dead_regions_are_not_active(geo):
    """The boundary between two sensors and the outer border must be inactive."""
    _, gm, lay, z = geo
    y = lay.globl(0)

    # midpoint between every pair of adjacent wafer columns
    for col in range(lay.cols - 1):
        x = (lay.wafer_centre(col) + lay.wafer_centre(col + 1)) / 2.0
        assert not _is_active(gm, x, y, z), f"boundary after column {col} is active"

    # outer border, just inside the plane edge
    edge = lay.dim / 2.0 - 0.005
    assert not _is_active(gm, -edge, y, z), "outer border is active"


def test_the_dead_cross_is_exactly_two_sensor_margins(geo):
    """Where the dead region starts and stops, to the micron.

    The rim is the whole of the dead region between two butted sensors, so a
    wrong margin does not show up as a hole in the pad grid -- the pads all stay
    where they belong -- but as a dead cross of the wrong width. Scanning across
    the boundary is what pins it down.
    """
    _, gm, lay, z = geo
    y = lay.globl(0)
    eps = 0.0002  # 2 um, in cm

    for col in range(lay.cols - 1):
        mid = (lay.wafer_centre(col) + lay.wafer_centre(col + 1)) / 2.0
        half_dead = lay.margin + lay.gap / 2.0
        assert not _is_active(gm, mid - half_dead + eps, y, z), (
            f"boundary {col}: still dead {half_dead * CM_TO_MM:.3f} mm before the middle")
        assert _is_active(gm, mid - half_dead - eps, y, z), (
            f"boundary {col}: dead region reaches past the rim on the low side")
        assert _is_active(gm, mid + half_dead + eps, y, z), (
            f"boundary {col}: dead region reaches past the rim on the high side")


def test_the_sensitive_area_is_the_sensor_minus_its_rim(geo):
    """The sensitive volume must be the pad array, not the whole sensor.

    If it were the whole sensor, a hit in the rim would be assigned to a pad --
    and, worse, to the *neighbouring* sensor's index range, because the
    segmentation grid does not stop at the edge of the volume it divides.
    """
    _, _, lay, _ = geo
    assert abs(lay.active - (lay.wafer - 2.0 * lay.margin)) < 1e-9
    assert abs(lay.active - lay.pads_per_wafer * lay.pad) < 1e-9
    # The pixel pitch is the diode plus the inter-pad margin. That margin is
    # silicon belonging to no diode and is deliberately not a volume (it is
    # 0.36% of the area and far below any step length), so it is checked here on
    # the constants rather than by navigating the geometry.
    assert abs(lay.pad - (lay.diode + lay.pad_margin)) < 1e-9
    assert lay.diode < lay.pad, "the active pad must be smaller than the pixel pitch"


def test_cellid_round_trip_gives_global_indices(geo):
    """position -> cellID -> position, with x/y running 0..npads-1 over the plane."""
    import ROOT

    det, _, lay, _ = geo
    seg = det.readout(READOUT).segmentation()
    dec = det.readout(READOUT).idSpec().decoder()
    shift = next(s for s in range(64) if dec.get(1 << s, "wafer") == 1)

    ids, bad = set(), []
    for i in range(lay.npads):
        for j in range(lay.npads):
            vol_id = lay.wafer_id(i, j) << shift
            local = ROOT.dd4hep.Position(lay.local(i), lay.local(j), 0)
            world = ROOT.dd4hep.Position(lay.globl(i), lay.globl(j), 0)
            cid = seg.cellID(local, world, vol_id)
            ids.add(cid)
            if (dec.get(cid, "x"), dec.get(cid, "y")) != (i, j):
                bad.append((i, j, dec.get(cid, "x"), dec.get(cid, "y")))
            back = seg.position(cid)
            if abs(back.X() - lay.local(i)) > 1e-9 or abs(back.Y() - lay.local(j)) > 1e-9:
                bad.append((i, j, back.X(), back.Y()))

    assert not bad, f"{len(bad)} cellID mismatches, e.g. {bad[:5]}"
    assert len(ids) == lay.npads ** 2, "cell IDs are not unique"


def test_no_overlaps(geo):
    """The wafers must not overlap each other nor stick out of the plane."""
    _, gm, _, _ = geo
    gm.CheckOverlaps(0.0001)
    overlaps = [o.GetName() for o in gm.GetListOfOverlaps()]
    assert not overlaps, f"{len(overlaps)} overlaps, e.g. {overlaps[:3]}"


def test_pad_positions_match_the_real_pad_map(geo):
    """The simulated pads must sit exactly where the test-beam pad map says.

    This is what catches a wrong sensor margin: ChannelMapper matches hits to
    channels by position, so a systematic offset is silently absorbed by its
    nearest-neighbour search instead of failing loudly.

    Agreement has to be exact here, to the rounding of the map file: the map is
    generated from these same constants, so anything else means the two have
    drifted apart and the map needs regenerating.
    """
    _, _, lay, _ = geo
    if lay.asus != 1:
        pytest.skip("pad map only covers a single ASU")
    if not os.path.exists(PAD_MAP):
        pytest.skip(f"pad map not found: {PAD_MAP}")

    columns = set()
    for n, line in enumerate(open(PAD_MAP)):
        parts = line.split()
        if n == 0 or len(parts) < 6:
            continue
        try:
            columns.add(round(float(parts[4]), 4))
        except ValueError:
            continue
    real = sorted(columns)

    assert len(real) == lay.npads, f"pad map has {len(real)} columns, geometry has {lay.npads}"
    worst = max(abs(lay.globl(i) * CM_TO_MM - real[i]) for i in range(lay.npads))
    assert worst < PAD_MAP_TOL_MM, f"pads are off the real map by up to {worst:.4f} mm"
    # The two tables must also agree end to end, or the wafer array as a whole
    # would be the wrong size and the disagreement would not be a rounding of
    # the pitch but a real geometry error.
    span_geo = (lay.globl(lay.npads - 1) - lay.globl(0)) * CM_TO_MM
    span_map = real[-1] - real[0]
    assert abs(span_geo - span_map) < 0.1, (
        f"the pad array spans {span_geo:.3f} mm, the map {span_map:.3f} mm")

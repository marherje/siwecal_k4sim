"""
Tests that the two visualisers draw the segmented silicon where it really is.

``analysis/tests/test_wafer_geometry.py`` proves geometry, readout and pad map
agree. That leaves one gap: the displays could still paint a continuous plane
over the guard rings, which is exactly the kind of error that goes unnoticed --
the picture looks fine, and hits that fall next to a dead region look like they
sit on active silicon.

So both visualisers are checked against the same source of truth, the compact
file:

  * ``event_display``  reads the wafer placements out of DD4hep directly
    (``use_geometry_placement``), so the test asserts it finds the wafers rather
    than the enclosing plane.
  * ``event_viewer``   has no access to the compact file at run time -- it works
    from the pad map -- so the test asserts the wafer rectangles it *infers*
    match the DD4hep wafers to well under a pad pitch.

The units are the usual trap: DD4hep works in cm, the pad map and the viewer in
mm. Everything is compared in mm.

Run with:
    source init_key4hep.sh   # release pinned in .key4hep-release
    export LD_LIBRARY_PATH=$PWD/install/lib64:$PWD/install/lib:$LD_LIBRARY_PATH
    python -m pytest analysis/tests/test_viewer_wafer_display.py -v
"""

from __future__ import annotations

import json
import os
import sys

import pytest

REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
COMPACT = os.path.join(REPO_ROOT, "simulation", "geometry", "SND_compact.xml")
DISPLAY_DIR = os.path.join(REPO_ROOT, "event_display")
DISPLAY_CONFIG = os.path.join(DISPLAY_DIR, "detector_config.json")
MAPPINGS = os.path.join(REPO_ROOT, "mappings")

CM_TO_MM = 10.0

# The viewer places an active-area edge half a pad pitch outside the outermost
# pad centre, which is exact because the pads butt against each other inside a
# sensor. Allow a tenth of a millimetre for the rounding of the map file.
EDGE_TOL_MM = 0.1

# Inactive rim of a sensor [mm], Ecal_WaferMarginX. Repeated here so that the
# expected size of the dead cross is written out rather than inferred.
WAFER_MARGIN_MM = 0.61


def _require():
    if not os.path.exists(COMPACT):
        return f"compact file not found: {COMPACT}"
    try:
        import ROOT  # noqa: F401
    except ImportError as exc:
        return f"key4hep not set up: {exc}"
    try:
        import plotly  # noqa: F401
    except ImportError as exc:
        return f"viewer venv not active: {exc}"
    return None


pytestmark = pytest.mark.skipif(_require() is not None, reason=_require() or "")


@pytest.fixture(scope="module")
def display():
    """``(boxes_cm, dd4hep detector)`` for the silicon entry of the display.

    Both come from the same single ``fromXML`` call: ``dd4hep.Detector`` is a
    process-wide singleton, and loading the compact file into it twice hangs.
    """
    sys.path.insert(0, DISPLAY_DIR)
    from event_display_eve import extract_z_from_geometry
    import ROOT

    with open(DISPLAY_CONFIG) as handle:
        config = json.load(handle)
    extract_z_from_geometry(COMPACT, config)  # loads the geometry
    entry = next(g for g in config["geometry"] if g["geo_extract"].get(
        "use_geometry_placement") and "Si" in g["name"])
    boxes = entry.get("boxes_cm")
    assert boxes, "the silicon entry extracted no placed volumes"
    return boxes, ROOT.dd4hep.Detector.getInstance()


@pytest.fixture(scope="module")
def display_boxes(display):
    return display[0]


@pytest.fixture(scope="module")
def viewer_model():
    """A ``DetectorModel`` built from the production pad maps."""
    sys.path.insert(0, REPO_ROOT)
    from event_viewer._geometry import DetectorGeometry
    from event_viewer._pad_map import PadMap
    from event_viewer.model.detector import DetectorModel

    pad_map = PadMap.from_files(
        {"default": "fev10_rotate_chip_channel_x_y_mapping.txt",
         12: "fev11_cob_good_rotate_chip_channel_x_y_mapping.txt"},
        base_dir=MAPPINGS)
    geometry = DetectorGeometry()
    return DetectorModel(geometry, pad_map, geometry.slab_z_mm,
                         geometry.slab_w_thickness_mm)


def _dd4hep_wafers_mm(boxes, margin_mm=0.0):
    """``(x0, x1, y0, y1)`` in mm of the distinct sensors, from one layer.

    The display draws the physical sensor; pass ``margin_mm`` to shrink the
    rectangles to the active area inside it, which is what the viewer draws.
    """
    z0 = min(b[2] for b in boxes)
    return sorted(
        ((x - dx) * CM_TO_MM + margin_mm, (x + dx) * CM_TO_MM - margin_mm,
         (y - dy) * CM_TO_MM + margin_mm, (y + dy) * CM_TO_MM - margin_mm)
        for x, y, z, dx, dy, _dz in boxes if z == z0)


# ------------------------------------------------------------ event_display --
def test_display_draws_wafers_not_planes(display):
    """The silicon boxes are wafer-sized, not plane-sized.

    A regression here means the display fell back to one box per layer, which is
    what hid the dead regions before the geometry was segmented.
    """
    display_boxes, det = display
    wafer_half = det.constantAsDouble("Ecal_WaferSizeX") / 2.0
    plane_half = det.constantAsDouble("Ecal_dim_x") / 2.0

    halves = {round(b[3], 6) for b in display_boxes}
    assert halves == {round(wafer_half, 6)}, (
        f"silicon half-sizes {halves} cm are not the wafer's {wafer_half} cm "
        f"(the plane would be {plane_half} cm)")


def test_display_covers_every_wafer_of_every_layer(display):
    """One box per wafer per layer, all layers, no duplicates."""
    display_boxes, det = display
    n_wafers = (int(det.constantAsDouble("Ecal_NASUsX"))
                * int(det.constantAsDouble("Ecal_NWafersX"))
                * int(det.constantAsDouble("Ecal_NASUsY"))
                * int(det.constantAsDouble("Ecal_NWafersY")))
    n_layers = int(det.constantAsDouble("Ecal_NLayers"))

    assert len(display_boxes) == len(set(display_boxes)), "duplicate boxes drawn"
    assert len(display_boxes) == n_wafers * n_layers
    per_layer = {}
    for box in display_boxes:
        per_layer.setdefault(box[2], []).append(box)
    assert len(per_layer) == n_layers
    assert all(len(v) == n_wafers for v in per_layer.values())


# ------------------------------------------------------------- event_viewer --
def test_viewer_infers_the_wafers_from_the_pad_map(viewer_model, display_boxes):
    """The rectangles drawn by the viewer are the DD4hep wafers, in mm.

    The viewer never reads the compact file: it recovers the wafers from the
    steps in the pad map, then adds the sensor rim. This is the check that the
    inference is right, and it is what would catch a mapping file whose gap no
    longer matches the simulated geometry. See EDGE_TOL_MM for why the two are
    compared to half a millimetre rather than exactly.
    """
    expected = _dd4hep_wafers_mm(display_boxes, margin_mm=WAFER_MARGIN_MM)
    got = sorted(viewer_model.wafer_rects)

    assert len(got) == len(expected), (
        f"viewer draws {len(got)} wafers, geometry has {len(expected)}")
    for rect, ref in zip(got, expected):
        for value, reference, edge in zip(rect, ref, ("x0", "x1", "y0", "y1")):
            assert abs(value - reference) < EDGE_TOL_MM, (
                f"wafer edge {edge}: viewer {value:.3f} mm vs "
                f"geometry {reference:.3f} mm")


def test_viewer_leaves_the_boundary_between_sensors_empty(viewer_model):
    """No drawn rectangle covers the dead cross, and the cross is the right size.

    The sensors are butted, so the whole of the dead region between two of them
    is their two 0.61 mm rims: 1.22 mm, and the pads either side are 6.75 mm
    apart against the 5.53 mm pitch inside a sensor.
    """
    rects = sorted(viewer_model.wafer_rects)
    xs = sorted({(r[0], r[1]) for r in rects})
    assert len(xs) > 1, "expected more than one wafer column"

    gap = xs[1][0] - xs[0][1]
    assert gap > 0, "sensors overlap; the dead cross would be painted over"
    expected = 2 * WAFER_MARGIN_MM
    assert abs(gap - expected) < EDGE_TOL_MM, (
        f"dead cross is {gap:.3f} mm, expected {expected:.3f}")


def test_pad_cells_are_one_pitch_and_never_overlap(viewer_model):
    """Each pad occupies exactly one pitch-wide cell, the guard ring its own.

    This is what a scatter marker cannot give: its size is in pixels, so at the
    4x4 grid's ~150 px per layer the 9 px hit markers were more than twice the
    4 px a 5.5 mm pad spans, overlapping their neighbours and reaching over the
    guard ring. Cells are in data coordinates and stay exact at any zoom.
    """
    import numpy as np

    pitch = viewer_model.pad_pitch
    for slab in (0, 12):  # 12 uses the other mapping file
        x_edges, y_edges, x_index, y_index = viewer_model.pad_cell_grid(slab)
        pads = viewer_model.pads_for_slab(slab)

        for edges, name in ((x_edges, "x"), (y_edges, "y")):
            widths = np.diff(edges)
            assert (widths > 0).all(), f"{name} cell edges are not increasing"
            pad_cells = np.isclose(widths, pitch)
            # Everything that is not a pad is a dead region, and there must be
            # at least one of those or the sensors were not separated at all.
            gaps = widths[~pad_cells]
            assert gaps.size >= 1, f"no dead-region cell along {name}"
            assert np.allclose(gaps, 2 * WAFER_MARGIN_MM, atol=1e-3), (
                f"{name} dead-region cells are {gaps} mm wide")
            assert pad_cells.sum() == len(np.unique(np.round(
                pads[:, 0 if name == "x" else 1], 3)))

        # Every pad lands on a cell of its own: no two pads share one.
        cells = {(y_index[round(float(y), 3)], x_index[round(float(x), 3)])
                 for x, y in pads}
        assert len(cells) == len(pads), "two pads were placed on the same cell"


def test_scene3d_draws_one_trace_per_plane(viewer_model):
    """Wafers are batched per plane, so the scene keeps its original trace count.

    Splitting a plane into wafers must not multiply the number of plotly traces:
    the accumulated-cluster scenes are already heavy, and a 4x trace count there
    is a visible slowdown in the browser.
    """
    from event_viewer.viz.scene3d import DetectorScene3D

    n_layers = viewer_model.slab_z_mm.size
    n_w = len(viewer_model.tungsten_boxes())
    traces = DetectorScene3D(viewer_model).base_figure().data
    assert len(traces) == n_layers + n_w

    # ...and each silicon trace really does carry all the wafers.
    n_wafers = len(viewer_model.wafer_rects)
    assert len(traces[0].x) == 4 * n_wafers

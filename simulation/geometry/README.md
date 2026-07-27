# Detector geometry

DD4hep compact description of the SiW-ECAL test beam setup, plus the two scripts
that generate and verify it.

| file | what it is |
|---|---|
| `SND_compact.xml` | top-level compact file — constants, readout, world/beam setup. This is the file every job points at. |
| `SiPadDetector.xml` | the detector element itself: layer stack, wafer tiling, sensitive volumes. Built by the `SiPadDetector` plugin in `detector_plugin/`. |
| `elements.xml`, `materials.xml` | GDML element table and the custom materials (glue, PCB, …). |
| `make_readout.py` | generates the `<readout name="SiPadHits">` block for a given number of ASUs. |
| `dead_zone_test.py` | fires muons through the inter-wafer gaps to check Geant4 really drops the energy there. |
| `general_tests/` | scaled-up 36 cm and 54 cm variants, kept as worked examples — see its own README. |

The production geometry is the single-ASU test-beam layout: 15 layers, 2x2
wafers per layer, 16x16 pads per wafer, so 32x32 pads and 18x18 cm per layer.

---

## Segmentation in one paragraph

The readout uses a `MultiSegmentation` keyed on the `wafer` physVolID, with one
`TiledLayerGridXY` entry per wafer. Each entry has the same grid (the pixel pitch
`Ecal_CellSizeX` = pad + inter-pad margin) but a different offset, shifted by a
whole pad array, so that the decoded `x`/`y` come out as **global** pad indices
over the layer (0..31 here) rather than per-wafer indices. That is why the block
has to be generated: 1 ASU needs 4 entries, 3x3 ASUs need 36.

The offsets are expressed in the frame of the sensitive volume, so the shift is
the active pad array `Ecal_NPadsPerWaferX*Ecal_CellSizeX`, **not** the physical
sensor size — the 0.61 mm inactive rim plays no part in the segmentation.

---

## `make_readout.py`

```bash
python simulation/geometry/make_readout.py            # 1x1 ASU (default)
python simulation/geometry/make_readout.py 2 2        # 2x2 ASUs
python simulation/geometry/make_readout.py 3 3        # 3x3 ASUs
```

Prints the whole `<readout name="SiPadHits">` block to stdout. To change the
detector size:

1. set `Ecal_NASUsX` / `Ecal_NASUsY` in `SND_compact.xml` to the same values;
2. paste the script output over the existing `<readout>` block;
3. re-run the checks below.

`WAFERS_PER_ASU` and `PADS_PER_WAFER` at the top of the script must stay in sync
with `Ecal_NWafersX/Y` and `Ecal_NPadsPerWaferX/Y`. The `wafer` field is last in
the `<id>` encoding precisely so that widening it never moves the bit positions
of `system`/`layer`/`slice`/`x`/`y`.

---

## `dead_zone_test.py`

The static checks compare geometry against readout; they cannot show that Geant4
actually discards energy deposited in a gap — that depends on which volume
carries the sensitive detector, which is the easy thing to get wrong.

The script fires 10 GeV muons along z at a fixed (x, y) and measures the *plane
efficiency*: the fraction of (event, layer) pairs with at least one hit. A MIP
crossing active silicon fires every plane, so a pad centre must give 100% and a
dead region must collapse. Needs a working `ddsim`, so source the environment
first:

```bash
source init_siwecal_soft.sh
export LD_LIBRARY_PATH=$PWD/install/lib64:$PWD/install/lib:$LD_LIBRARY_PATH

python simulation/geometry/dead_zone_test.py                     # default scan
python simulation/geometry/dead_zone_test.py --x 0 --y -58.8     # single point
python simulation/geometry/dead_zone_test.py \
       --compact simulation/geometry/general_tests/SND_compact_36cm.xml
```

| flag | default | meaning |
|---|---|---|
| `--compact` | `SND_compact.xml` | geometry to test |
| `--events` | 20 | muons per target point |
| `--layers` | 15 | planes expected per event, the efficiency denominator |
| `--x`, `--y` | — | test one point instead of the default scan |

The default scan walks from a pad centre in to the middle of the dead cross. The
efficiency has to collapse between 0.7 mm and 0.5 mm off the middle, where the
two 0.61 mm rims start. The 0.7 mm row reads as a transition (~80%), not a
failure: the aim point is on pad 15 but only 0.09 mm from the rim, and a muon
scatters further than that before reaching the back of the stack. The
pad-boundary row is the control — inside a sensor pads butt against each other,
so a muon on a pad boundary still fires every plane, it just shares the charge.
Without it you cannot tell a correct dead region from a readout that silently
drops hits at every pad boundary.

Output goes to a fresh `mkdtemp` directory (path printed at startup) holding one
ROOT file and one `ddsim` log per target point; nothing is cleaned up, so delete
it when done.

---

## Checking a geometry change

Static — geometry vs readout vs pad map (the pad-map test skips itself for more
than one ASU, since `mappings/` only covers 32x32):

```bash
python -m pytest analysis/tests/test_wafer_geometry.py -v
```

Dynamic — the muon scan above. Between them they cover overlaps, pad centres
landing in silicon, cell ID uniqueness, and hits actually vanishing in the gaps.

Note that a compact file living in a subdirectory (as in `general_tests/`) needs
`../` in front of its `<gdmlFile>` and `<include>` refs.

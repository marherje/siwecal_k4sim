# Alternative detector sizes

Scaled-up variants of `../SND_compact.xml`, kept here as worked examples of the
two-level tiling (ASUs per layer x wafers per ASU). They are **not** used by the
production pipelines, which run the single-ASU test-beam geometry.

| file | ASUs | wafers | pads/layer | layer size | `wafer` field | MultiSeg entries |
|---|---|---|---|---|---|---|
| `../SND_compact.xml` | 1x1 | 2x2 | 32x32 | 18x18 cm | 2 bits | 4 |
| `SND_compact_36cm.xml` | 2x2 | 4x4 | 64x64 | 36x36 cm | 4 bits | 16 |
| `SND_compact_54cm.xml` | 3x3 | 6x6 | 96x96 | 54x54 cm | 6 bits | 36 |

In all three, `x`/`y` stay **global** pad indices over the whole layer (0..31,
0..63, 0..95) and keep the exact same bit positions — only `wafer`, which sits
last in the encoding, gets wider.

The dead regions are of two kinds and they are not the same width:

- **2.1 mm** between wafers inside an ASU (the guard ring, measured from
  `mappings/fev10_rotate_chip_channel_x_y_mapping.txt`);
- **1.9 mm** between wafers of adjacent ASUs, which is not a free parameter: the
  wafers span 178.1 mm inside a 180 mm ASU, so butting ASUs at `Ecal_ASUPitch`
  leaves 2 x 0.95 mm. Raise the pitch if the real modules are mounted with extra
  mechanical clearance.

## Making another size

Set `Ecal_NASUsX` / `Ecal_NASUsY`, then regenerate the readout block, which needs
one entry per wafer because each wafer carries its own segmentation offset:

```bash
python simulation/geometry/make_readout.py 2 2     # 36 cm
python simulation/geometry/make_readout.py 3 3     # 54 cm
```

Paste the output over the `<readout name="SiPadHits">` block. Remember that a
compact file in this folder needs `../` in front of its `<gdmlFile>` and
`<include>` refs.

## Checking a size

Static checks — geometry vs readout vs pad map (the pad-map test skips itself
when there is more than one ASU, since the map only covers 32x32):

```bash
python -m pytest analysis/tests/test_wafer_geometry.py -v
```

Dynamic check — does Geant4 actually drop hits in the dead regions:

```bash
python simulation/geometry/dead_zone_test.py \
       --compact simulation/geometry/general_tests/SND_compact_36cm.xml
```

Both were run against these two files: 0 overlaps, all pad centres in silicon,
cell IDs unique, and a muon aimed at global pad (60,60) of the 36 cm layer — the
far corner ASU — gives 300/300 hits in cell (60,60).

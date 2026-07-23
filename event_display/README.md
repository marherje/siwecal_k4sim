# SND Event Display

Interactive 3D event display for SND@LHC based on ROOT TEve.
Detector geometry and hit-display parameters are driven entirely by
`detector_config.json` — no positions are hardcoded in the Python script.

## Prerequisites

- key4hep environment sourced
- Simulation pipeline completed: `ShipHits.root` must exist in the working
  directory (produced by `job5_rntuple.py`)
- DD4hep geometry built: `../simulation/geometry/SND_compact.xml` must be present

---

## How to run the event display

```bash
python event_display_eve.py \
    --hits  ShipHits.root \
    --config detector_config.json \
    --window 0
```

| Option | Description |
|--------|-------------|
| `--hits` | Path to `ShipHits.root` (RNTuple output of `job5_rntuple.py`) |
| `--config` | Path to the JSON config (default: `detector_config.json`) |
| `--window` | Time-window / event index to display (default: `0`) |

Close the TEve window to exit.

---

## Config file structure (`detector_config.json`)

```json
{
  "detectors": [
    {
      "name":        "SiTarget_StripX",
      "ntuple":      "SiTargetMeas",
      "filter":      {"plane": 0},
      "color":       [1.0, 0.4, 0.7],
      "voxel":       {"x": 0.003775, "y": 4.9, "z": 0.015},
      "layers_z_cm": [ ... ]
    },
    ...
  ],
  "geometry": [
    {
      "name":         "SiTarget planes",
      "color":        [0.2, 0.4, 0.9],
      "transparency": 90,
      "voxel":        {"x": 20.0, "y": 20.0, "z": 0.015},
      "layers_z_cm":  [ ... ]
    },
    ...
  ]
}
```

- **`detectors`** — each entry reads one RNTuple from `ShipHits.root` and draws
  one box per hit; `filter` applies an equality cut on any integer field
  (e.g. `{"plane": 0}` keeps only StripX hits).
- **`geometry`** — transparent outline planes drawn in the global scene as a
  detector reference frame.
- **`voxel`** — half-sizes of the box drawn per hit, in cm.
- **`color`** — RGB triplet `[r, g, b]` in the range `[0, 1]`.
- **`geo_extract`** — fills the entry's layer positions from the DD4hep compact
  file instead of hardcoding them. `path_contains` selects the subdetector and
  `slice_material` (or `slice_index`) the slices to pick up. Two flags refine
  what is drawn:
  - `use_geometry_thickness` — take each layer's `dz` from the geometry rather
    than from `voxel.z`.
  - `use_geometry_placement` — draw one box per *placed volume*, at its own
    `(x, y, z)` and with its own half-sizes, ignoring `voxel` entirely. This is
    what the silicon uses: a sensitive slice is an array of 89.7 mm silicon
    sensors, each an inactive 0.61 mm rim around 16x16 pads of 5.53 mm pitch, so
    a single full-plane box per layer would hide the dead regions. When a rule
    matches twice along the same path -- the sensor and the pad array inside it
    are both silicon -- only the outermost volume is drawn, i.e. the sensor.

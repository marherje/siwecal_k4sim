# DD4hep Detector Plugins

## Overview

Three sub-detectors, each a C++ DD4hep plugin (`dd4hep_add_plugin` in CMakeLists.txt):

| Plugin target | Source | Readout | Segmentation |
|--------------|--------|---------|--------------|
| `SiTargetDetector_plugin` | `detector_plugin/SiTargetDetector.cpp` | `SiTargetHits` | MultiSegmentation(CartesianStripX/Y keyed on `plane`) |
| `SiPadDetector_plugin` | `detector_plugin/SiPadDetector.cpp` | `SiPadHits` | CartesianGridXY |
| `MTCDetector_plugin` | `detector_plugin/MTCDetector.cpp` + `CartesianStripXStereo.cpp` | `MTCDetHits` | MultiSegmentation(CartesianStripXStereo U/V + CartesianGridXY keyed on `plane`) |

---

## Compact XML Structure

Top-level: `simulation/geometry/SND_compact.xml`  
Includes sub-detectors via `<include ref="..."/>`.

Detector ordering (Z along beam axis, particle enters at Z ≈ −1000 mm):
1. `SiTarget.xml` — upstream silicon strip tracker (20 layers, X/Y strips)
2. `SiPadDetector.xml` — silicon pad ECAL (middle, 20 layers)
3. `MTCDetector.xml` — downstream muon/tau catcher (3 stations × 15 layers)

### Readout declaration pattern
```xml
<readout name="MTCDetHits">
  <segmentation type="MultiSegmentation" key="plane">
    <segmentation type="CartesianStripXStereo" key_value="0"
                  strip_size_x="MTC_SciFiChannelSize"
                  offset_x="-(MTC_max_env_width/2.0)+0.001*mm"
                  stereo_angle="5.0"/>
    <segmentation type="CartesianStripXStereo" key_value="1"
                  strip_size_x="MTC_SciFiChannelSize"
                  offset_x="-(MTC_max_env_width/2.0)+0.001*mm"
                  stereo_angle="-5.0"/>
    <segmentation type="CartesianGridXY" key_value="2"
                  grid_size_x="MTC_ScintChannelSize" grid_size_y="MTC_ScintChannelSize"/>
  </segmentation>
  <id>system:8,station:2,layer:8,slice:4,plane:2,strip:14,x:9,y:9</id>
</readout>
```

**Critical:** Custom segmentation attributes (`stereo_angle`, etc.) must be **literal numbers** — DD4hep does NOT evaluate constant expressions for custom plugins (uses raw `stringstream` parsing).

---

## SiTarget Detector (`SiTargetDetector.cpp`)

- 20 layers; each layer has two sensitive planes
- `plane=0` → CartesianStripX (X-measuring strips)
- `plane=1` → CartesianStripY (Y-measuring strips)
- MultiSegmentation keyed on `plane` bit (offset = 8+2+8+4 = 22 for `system:8,station:2,layer:8,slice:4,plane:2,...`)

---

## SiPad Detector (`SiPadDetector.cpp`)

- 20 layers, CartesianGridXY segmentation
- Used as ECAL; also feeds ACTS tracking

---

## MTC Detector (`MTCDetector.cpp`)

Three stations: MTC40 (half-height 200mm), MTC50 (250mm), MTC60 (300mm).

Each station: 15 layers. Each layer:

| Slice | Material | `plane` | Sensitive |
|-------|----------|---------|-----------|
| Iron absorber ×2 | Fe | — | no |
| SciFi U | Silicon | 0 | yes |
| Air gap | — | — | no |
| SciFi V | Silicon | 1 | yes |
| Scintillator | Scintillator | 2 | yes |

All sensitive slices share one readout (`MTCDetHits`) routed by `plane` field.

### CellID layout
`system:8, station:2, layer:8, slice:4, plane:2, strip:14, x:9, y:9`

Bit offsets:
- system: 0–7
- station: 8–9
- layer: 10–17
- slice: 18–21
- **plane: 22–23** ← routing key for MultiSegmentation and EventShuffler

---

## CartesianStripXStereo Segmentation (`detector_plugin/CartesianStripXStereo.cpp`)

Custom DD4hep segmentation for stereo SciFi strips.

**Strip assignment:**
```
strip = floor((local_x - local_y * tan(stereo_angle_deg * π/180)) / pitch)
```
- U plane (stereo_angle = +5°): strip shifts right as y increases
- V plane (stereo_angle = −5°): strip shifts left as y increases

**`position(cellID)`** returns strip centre at y=0 (x₀ = strip_centre_at_y=0, y=0, z=0). This gives the reference stereo coordinate; the actual y of a hit is stored separately in `hit.position.y` by `SND_SciFiAction`.

**Implementation class naming:** Plugin name (`CartesianStripXStereo`) must differ from impl class name (`CartesianStripXStereoImpl`) because `DECLARE_SEGMENTATION` creates a stub struct named after the plugin string.

**Unit:** `position()` returns ROOT cm (not mm). See unit convention table in the `key4hep-snd` skill.

---

## DDG4 Sensitive Action (`simulation/ddg4/SND_SciFiAction.cpp`)

Replaces the default DDG4 calorimeter action for the MTC detector.

Per Geant4 step:
1. `process()` — computes stereo-corrected strip cellID via `m_segmentation.cellID(step)`. Sets `hit.position.x = seg_p.x()` (strip centre at y=0, ROOT cm). Accumulates energy-weighted step-Y: `sumEY[cellID] += edep * y_mid`.
2. `end()` — overwrites `hit.position.y = sumEY / sumE` (energy-weighted avg Y in Geant4 mm).

Activated in ddsim steering:
```python
SIM.action.mapActions["MTC"] = "SND_SciFiAction"
```
String `"MTC"` must match the detector name in `SND_compact.xml`.

**Build requirement:** Links `DD4hep::DDG4` (transitive Geant4 deps) → must have `INSTALL_RPATH_USE_LINK_PATH TRUE`.

For full details on position encoding and unit inconsistencies, see `docs/mtc_scifi_hit_pipeline.md`.

---

## Adding a New Detector Plugin

1. Create `detector_plugin/MyDetector.cpp`
2. Add to CMakeLists.txt:
   ```cmake
   dd4hep_add_plugin(MyDetector_plugin SOURCES detector_plugin/MyDetector.cpp)
   target_link_libraries(MyDetector_plugin DD4hep::DDCore)
   ```
3. Add `<include ref="MyDetector.xml"/>` to `SND_compact.xml`
4. `/build` to rebuild
5. Update `init_key4ship.sh` if a new env variable is needed

# ACTS Track Reconstruction Integration

## Overview

ACTS (A Common Tracking Software) is integrated via `ACTSGeoSvc` (a Gaudi Service) and `ACTSProtoTracker` (a Gaudi Algorithm), both compiled into `libSND_reco.so`.

ACTS version: 44.3.0 (from key4hep 2026-02-01 CVMFS stack).

---

## ACTSGeoSvc

**Files:** `gaudi_source/ACTSGeoSvc.h`, `gaudi_source/ACTSGeoSvc.cpp`  
**Interface:** `gaudi_source/ISNDGeoSvc.h` (extends `IActsGeoSvc` from k4ActsTracking)

### Purpose
Loads DD4hep geometry, walks the TGeo tree, and builds ACTS `TrackingGeometry` consisting of `PlaneSurface` objects for each sensitive detector plane.

### Configuration
```python
from Configurables import ACTSGeoSvc
geo = ACTSGeoSvc("ACTSGeoSvc")
geo.CompactFile    = "simulation/geometry/SND_compact.xml"  # required
geo.MTCStereoAngle = 5.0  # degrees, must match CartesianStripXStereo XML value
ApplicationMgr(..., ExtSvc=[geo])  # MUST be in ExtSvc, not TopAlg
```

### Interface methods
```cpp
// In algorithm initialize():
ServiceHandle<ISNDGeoSvc> m_geoSvc{this, "ACTSGeoSvc", "ACTSGeoSvc"};

const Acts::TrackingGeometry& geo  = m_geoSvc->trackingGeometry();
const Acts::GeometryContext&  gctx = m_geoSvc->geometryContext();
const std::vector<const Acts::Surface*>& surfs = m_geoSvc->allSurfaces();
const Acts::Surface* surf = m_geoSvc->surfaceByAddress(detID, station, layer, plane);
```

`surfaceByAddress(detID, station, layer, plane)` is the primary lookup used by measurement converters.

### Geometry construction (ACTSGeoSvc.cpp)
1. Loads DD4hep geometry via `dd4hep::Detector::getInstance().fromXML(...)`
2. Walks `gGeoManager` TGeo tree, finds sensitive volumes by detector name and parses `(detID, station, layer, plane)` from the volume name
3. Extracts global Z position and half-sizes from TGeo matrices (TGeo in **cm** → multiply ×10 for ACTS mm)
4. Constructs `Acts::PlaneSurface` with `Acts::RectangleBounds`. Surface rotation is `rot90Y` for SiTarget/SiPad and `R_X(±α) · rot90Y` for MTC SciFi (`+α` on U planes, `−α` on V planes, with `α = MTCStereoAngle`) — matches `CartesianStripXStereo` segmentation
5. Wraps surfaces in `Acts::PlaneLayer` → `Acts::TrackingVolume` → `Acts::TrackingGeometry`
6. Populates two parallel structures:
   - `m_allSurfaces` — flat vector sorted by ACTS-X (beam axis)
   - `m_surfaceByAddressMap` — `(detID, station, layer, plane) → const Acts::Surface*` for fast address lookup

### Surface material
Surface material is **required** for the CKF multiple-scattering and energy-loss correctors to do anything — without it, each iron slab is treated as vacuum and the fit absorbs MS kicks into q/p (this caused the 3–10× momentum bias documented in `gaudi_jobs/muon_pipeline/MTC_OVERCURVATURE_FINDINGS.md`, resolved 2026-05-17).

After surface construction, `ACTSGeoSvc` queries `dd4hep::Detector::material(name)` for X₀, L₀, ρ and attaches `Acts::HomogeneousSurfaceMaterial` to every surface:

| Detector | Material attached |
|---|---|
| SiTarget plane=0 | TungstenDens1910 (3.5 mm) + Silicon (0.3 mm) |
| SiTarget plane=1 | Silicon (0.3 mm) |
| SiPad | TungstenDens1910 (3.5 mm) + Silicon (0.65 mm) |
| MTC plane=0 (U) | Iron (53 mm) + Scintillator (1.35 mm) |
| MTC plane=1 (V) | Scintillator (1.35 mm) |

The thin Cu/CF layers are intentionally skipped.

### Unit notes
- TGeo positions are in **cm**; all values multiplied by 10 before passing to ACTS
- ACTS uses **mm** and **GeV** throughout
- `Acts::UnitConstants::mm` = 1.0, `Acts::UnitConstants::GeV` = 1.0

---

## Measurement Converters

The three converters write `edm4hep::TrackerHit3D` collections that `ACTSProtoTracker` later reads. Each converter encodes the address fields that `ACTSGeoSvc::surfaceByAddress(detID, station, layer, plane)` needs:

| Converter | Output collection | `type` | `quality` | `cellID` decoding done by tracker |
|---|---|---|---|---|
| `SiTargetMeasConverter` | `SiTargetMeasurements` | 0 | strip plane (0=X, 1=Y) | layer (not yet — see Issue 3 outstanding work) |
| `SiPadMeasConverter` | `SiPadMeasurements` | 1 | layer (0…19) | — |
| `MTCSciFiMeasConverter` | `MTCSciFiMeasurements` | 2 | plane (0=U, 1=V) | `station`, `layer` from `cellID` via `BitFieldCoder` |

### SiTargetMeasConverter (`gaudi_source/SiTargetMeasConverter.cpp`)
Converts `SiTargetHitsWindowed` (EDM4HEP `SimTrackerHit`) → `TrackerHit3D`. `type=0` (SiTarget), `quality=plane` (0=StripX, 1=StripY). Strip coordinate is the global x or y (depending on plane) stored in `position`.

### SiPadMeasConverter (`gaudi_source/SiPadMeasConverter.cpp`)
Converts `SiPadHitsWindowed` (`SimCalorimeterHit`) → `TrackerHit3D`. `type=1`, `quality=layer` (0…19) — used by `ACTSProtoTracker` for `surfaceByAddress(1, -1, layer, -1)`. 2D pad coordinate is `(position.x, position.y)`; per-coord variance = pitch²/12 with `PixelSizeX = PixelSizeY = 5.5 mm` by default.

### MTCSciFiMeasConverter (`gaudi_source/MTCSciFiMeasConverter.cpp`)
Converts `MTCSciFiHitsWindowed` (plane 0 = U, plane 1 = V) → 1D `TrackerHit3D` on MTC SciFi surfaces. Scintillator hits (plane 2) are skipped. Strip coordinate (`cos α · pos.x − sin α · pos.y` on a U plane, mirrored for V) is stored in `position.x` [mm]; variance = pitch²/12. `type=2`, `quality=plane`.

**Status:** Fully implemented and active in `muon_pipeline`. Not wired in other pipelines (see [Pipeline differences](#pipeline-differences--mtc-and-b-field)).

### SNDMeasurement struct (internal to ACTSProtoTracker)
```cpp
struct SNDMeasurement {
  const Acts::Surface* surface;
  double localCoord;    // 1D strip coordinate (SiTarget/SciFi) or X (SiPad)
  double localCoord2;   // Y coordinate (SiPad only, is2D=true)
  double variance;      // measurement uncertainty²
  double variance2;     // Y uncertainty² (SiPad only)
  bool   is2D;          // true for SiPad (pad detector)
  int    detectorID;    // 0=SiTarget, 1=SiPad, 2=MTC
  int    plane;         // strip plane index (or -1 for SiPad)
  float  time;          // hit time [ns]
  float  eDep;          // energy deposit [GeV]
};
```

---

## ACTSProtoTracker

**File:** `gaudi_source/ACTSProtoTracker.cpp`

### Algorithm flow
1. `initialize()` — retrieves `ACTSGeoSvc`, builds the `BitFieldCoder` for MTC `cellID` decoding, configures the `IronSlabBField` provider
2. `execute()` — per event:
   a. **Assemble measurements** by reading the three input collections and looking up each hit's surface through `m_geoSvc->surfaceByAddress(...)` (see [Measurement Converters](#measurement-converters) for the encoding scheme)
   b. `findSeeds()` — 2D Hough transform over assembled measurements
   c. For each seed: run the `CombinatorialKalmanFilter` (CKF) → `Acts::TrackContainer`
   d. Write output tracks to the `ACTSTracks` EDM4HEP collection
3. `finalize()` — log statistics

### Hough seeding
2D Hough transform over (x, z) or (y, z) plane. Parameters:
- `HoughBinSize` — accumulator bin width
- `HoughHalfSize` — accumulator range
- `HoughMinVotes` — min hits per bin to form a seed
- `AutoSeed` — if true, automatically constructs seeds from Hough peaks
- `MaxSeeds` — cap on number of seeds per event
- `SeedMomentum`, `SeedCharge` — qOverP prior. `SeedCharge = -1` matches the mu⁻ test gun; a wrong sign makes the CKF fight the B-field (see `MTC_OVERCURVATURE_FINDINGS.md`)

### CKF setup

The B-field provider is chosen per event based on the `IronFieldRanges` property:

- **`IronFieldRanges` set** (muon_pipeline): uses `IronSlabBField`, a custom `Acts::MagneticFieldProvider` that returns a configurable `By` inside each registered rectangular slab (the MTC outer iron absorbers) and zero everywhere else.
- **`IronFieldRanges` empty** (other pipelines): falls back to `Acts::ConstantBField(BFieldX, BFieldY, BFieldZ)`, all defaulting to 0.

The propagator uses a `DirectNavigator` wrapped in `SNDFixedNavigator`, which injects `m_geoSvc->allSurfaces()` (sorted by ACTS-X) into the navigator state. This forces the CKF to visit every sensitive surface in order, regardless of `TrackingVolume` boundaries — necessary because the SND geometry is built as a single flat surface list and a conventional `Navigator` cannot route between the surfaces correctly.

```cpp
auto stepper    = Acts::EigenStepper<>(bField);
auto navigator  = SNDFixedNavigator{};   // wraps Acts::DirectNavigator
auto propagator = Acts::Propagator(std::move(stepper), std::move(navigator));
Acts::CombinatorialKalmanFilterOptions ckfOptions(
    gctx, mctx, cref(cctx), ckfExtensions, pOptions,
    /*multipleScattering=*/true,
    /*energyLoss=*/true);
```

The chi² cut for `MeasurementSelector` is set via `ChiSquareCut` (default 15).

### Track output
Tracks are stored in an `edm4hep::TrackCollection` named `ACTSTracks`. Seeds (Hough candidates before fitting) are also stored in `ShipHits.root` via job5.

### Debug: `DIAG-TS` per-surface dump
When the message threshold is `INFO`, `ACTSProtoTracker` prints one row per CKF track state with fields:
```
[DIAG-TS] evt=… seed=… z=… <region> M=… O=… H=… pred(loc0,loc1) sig(σ0,σ1) [meas(loc0,loc1) cs=N]
```
- `M`/`O`/`H` = MeasurementFlag / OutlierFlag / HoleFlag (each 0 or 1; outliers also have `M=1`)
- `pred` and `sig` are the predicted (not smoothed) state at the surface
- `meas` and `cs` (= candidate source-link count) appear when a measurement was attached (inlier or outlier)

This dump is the primary diagnostic for per-detector inlier/outlier accounting. Used heavily in the Issue 3 investigation below.

---

## Magnetic Field

The MTC contains iron absorber slabs that carry a magnetic field (for muon momentum measurement). `ACTSProtoTracker` models this with a custom `IronSlabBField` provider that returns `By` inside each iron slab volume and zero elsewhere.

The slab ranges and field strength are passed via `IronFieldRanges` — a flat list of `[xlo, xhi, ylo, yhi, zlo, zhi, by]` septets in ACTS coordinates (mm, T). In `muon_pipeline/job4_tracking.py` these are computed automatically from `parse_geometry.SNDGeometry` (which reads `SND_compact.xml`), so the field map stays in sync with the geometry.

The SiTarget and SiPad regions have no magnetic field; `IronSlabBField` returns zero there automatically.

Other pipelines (no `IronFieldRanges`) fall back to `Acts::ConstantBField(0, 0, 0)`.

---

## Pipeline differences — MTC and B-field

`muon_pipeline` is the complete reference pipeline. It is the only one that wires MTC and the iron B-field:

| Feature | `muon_pipeline` | other pipelines |
|---|---|---|
| `MTCSciFiMeasConverter` in `TopAlg` | Yes | No |
| `InputMTC` wired to `ACTSProtoTracker` | Yes (`MTCSciFiMeasurements`) | No (default collection absent) |
| MTC U×V crossings in Hough seeding | Yes | No |
| MTC 1D hits fed into Kalman fitter | Yes | No |
| `IronFieldRanges` (MTC iron B-field) | Yes, from `parse_geometry` | No |
| Effective B-field | `IronSlabBField` (non-zero By in iron) | `ConstantBField(0,0,0)` |

To bring another pipeline up to full feature parity, follow `gaudi_jobs/muon_pipeline/job4_tracking.py` as the reference.

---

## Investigation findings (from `docs/remaining_acts_issues.md`)

This section records what we have learned about each of the four open ACTS issues. Items are presented in the order of the original list. Companion material: `MTC_OVERCURVATURE_FINDINGS.md` (the 2026-05-17 material-attachment investigation) and `known_issues.md`.

### Issue 1 — truth-track jumps in `diagnose_mtc_curvature.py`

**Symptom.** The visualised truth trajectory has occasional unrealistic jumps in coordinates at some MTC planes.

**What we have determined**
- The script reconstructs a global `(gx, gy)` per MTC layer in two modes:
  - **Paired** — when both U and V stereo states at the same layer have measurements, the script applies the stereo inversion `gy = (loc0_V − loc0_U) / (2·tan α)` to recover global y; gx comes from the local-x average. This is the correct formula given the U/V stereo rotation `R_X(±α) · rot90Y` in `ACTSGeoSvc`.
  - **Unpaired** — when only one of U or V has a measurement at a given layer (a missing SciFi station), the script falls back to writing the simulated `(x, y)` directly. This fallback is the source of the "weird jumps": the unpaired global y is the truth-particle y, while the paired global y is the U/V-strip-inverted estimate that includes any small smearing.
- A second contributor is a **sign-of-y discrepancy** between the DIAG-TS dump and the script plots — the script flips y for plotting (legacy convention from `event_display_eve.py`). This is purely cosmetic for the truth track but it is part of the y-sign trace below.

**Disposition.** The fallback behaviour is by design: the script chooses to plot *something* on layers with incomplete information rather than skip. Leave as-is. Document in the script's header that paired-layer points are stereo-inverted while unpaired-layer points are raw truth.

### Issue 2 — CKF end-of-track outliers

**Symptom.** Reconstructed tracks contain points at the end of the trajectory that lie far from the otherwise good track, and the CKF appears not to remove them.

**What we have determined**
- These are not measurement outliers in the statistical sense — they are **unpaired SciFi states** at the final MTC station, where one of the U or V planes recorded no hit. The CKF still creates a track state at the surface; `hasCalibrated()` returns true for both `InlierFlag` and `OutlierFlag` states (only true `HoleFlag` states return false), so both flow through the writer to the visualisation.
- They do **not degrade the CKF fit** — outliers are flagged and weighted accordingly, and the trajectory at upstream layers is unaffected.
- The DIAG-TS dump distinguishes them clearly: outliers show `M=1 O=1 H=0` with `cs ≥ 1` and visibly large `meas − pred` residuals (e.g. residuals ~40 mm on a σ ~8 mm prediction).

**Disposition.** Keep the visualisation as-is per user direction. If a cleaner display is later wanted, filter on `OutlierFlag` (or equivalently on a chi² threshold per state) at the writer step rather than in the CKF.

### Y-sign convention — cross-file trace

The user requested codebase-wide consistency on the y-sign convention even though the CKF itself is sign-agnostic. The convention is **global y is positive upward** (DD4hep / Geant4 default); the trace passes through six files:

| File | What it does | Sign |
|---|---|---|
| `simulation/ddg4/SND_SciFiAction.cpp` | Stores `position.x = seg_p.x() * 10` (cm→mm); y unchanged | + |
| `detector_plugin/CartesianStripXStereo.cpp` | `x_stereo = localPosition.X − localPosition.Y · tan(α)` | + |
| `gaudi_source/ACTSGeoSvc.cpp` | Surface rotation `R_X(+α) · rot90Y` on U planes, `R_X(−α) · rot90Y` on V planes | + (sign baked into rotation) |
| `gaudi_source/MTCSciFiMeasConverter.cpp` | `localCoord = cos(α) · pos.x` (stereo-projected) | + |
| `gaudi_source/ACTSProtoTracker.cpp` | Comment near line 1468 acknowledges `loc0 = −eBoundLoc0_geometric` (writer-side sign flip on output) | mixed |
| `event_display/event_display_eve.py` | Plots `−y` (cosmetic flip for display) | − for plotting |

**Findings that need fixing in comments / code**
1. `event_display/event_display_eve.py:464-465` — comment says "positive tilt → V plane"; this is **wrong**, V planes use the negative tilt. Comment only, no code change needed.
2. `gaudi_source/ACTSProtoTracker.cpp` — the smoothed `loc1_V` flips sign on a small fraction of states. This is consistent with the surface rotation but is surprising; document or sign-normalise at the writer.

The sign convention works end-to-end for the CKF; the artefacts are confined to comments and display.

### Issue 3 — combined-detector measurement assignment

**Symptom.** With SiTarget + SiPad + MTC all enabled, MTC reconstructed trajectories under- and over-estimate the bend in z-x and show flipped trajectories in z-y. MTC-only tracking does not show this.

**Root cause: broken gap-heuristic surface split.** Earlier versions of `ACTSProtoTracker::initialize()` partitioned `allSurfaces[]` into SiTarget / SiPad / MTC groups by finding the two largest z-gaps in the sorted surface list. With 40 + 20 + 90 = 150 surfaces, the **inter-station MTC40→MTC50 and MTC50→MTC60 gaps are larger than the inter-detector SiTarget→SiPad and SiPad→MTC gaps**, so the heuristic returned `SiTarget=60, SiPad=60, MTC=30` instead of the correct `40/20/90`.

Per-detector consequences:
- **MTC** — fine. Assembly uses `surfaceByAddress(2, station, layer, plane)` and never reads the gap-derived counts.
- **SiTarget** — works *by accident*. Nearest-z search runs over `allSurfaces[0..59]`, which contains 40 real SiTarget surfaces (z = −370…−155 mm) plus 20 SiPad surfaces (z = −140…+140 mm). Every SiTarget hit's z falls inside the SiTarget range, so the right surface always wins.
- **SiPad** — silently broken. Nearest-z search ran over `allSurfaces[60..119]`, which contained 60 MTC surfaces at z = 2200–3600 mm. Every SiPad hit (z ∈ [−140, +140] mm) was attached to an MTC surface ~2400 mm away. The CKF chi² then exceeded the `MeasurementSelector` cut and the measurement was silently dropped → **all 20 SiPad surfaces became holes**. Loss of the upstream y/x anchoring forced the MTC fit to absorb a bias.

**Diagnosis (DIAG-TS per-surface tally for `evt=0 seed=0`)**

| Region | n | Inliers (pre-fix) | Inliers (post-fix) | Holes (pre-fix) | Holes (post-fix) |
|--------|---|---|---|---|---|
| SiTarget | 40 | 40 | 40 | 0 | 0 |
| **SiPad** | 20 | **0** | **20** | **20** | **0** |
| MTC | 90 | 32 | 26 | 30 | 30 |

`nMeas` total: 72 → 86. Sum of CKF chi²: 66.3 → 144.9. The rise is dominated by the +14 SiPad inliers: with `PixelSize = 5.5 mm` the per-coord measurement σ is only 5.5/√12 ≈ 1.6 mm, and where the truth track crosses a pad boundary the pad-center residual reaches ~5 mm (≈ 3σ). chi²/n is **not** strictly worse — the two runs cover different measurement sets so the totals are not directly comparable.

**Fix applied**
- `gaudi_source/SiPadMeasConverter.cpp:199` — write the SiPad layer index into `TrackerHit3D::quality` (was hard-coded to `-1`).
- `gaudi_source/ACTSProtoTracker.cpp` SiPad branch (~lines 963–991) — replace the SiPad nearest-z search with `m_geoSvc->surfaceByAddress(1, -1, layer, -1)`.

**Outstanding work**
1. **SiTarget assembly still uses the nearest-z search bounded by `m_nSiTargetSurfaces`** (the gap-heuristic value, still 60 in the all-detector case). It works only by accident. Migrate to `surfaceByAddress(0, -1, layer, plane)` — requires writing the SiTarget layer into the converter output (re-purpose `quality` via packing, or decode `cellID` in the tracker).
2. **Once SiTarget is migrated, delete the gap-heuristic** and the `m_nSiTargetSurfaces` / `m_nSiPadSurfaces` members. The DIAG-TS region label currently uses them; switch to deriving the region from the surface's stored address.

**Still-pending verification.** The chi² jump 66 → 145 is consistent with the fix doing its job; what is **not yet verified** is whether the smoothed MTC trajectory matches truth better than before. Re-run `diagnose_mtc_curvature.py` and compare MTC |Δx| and `R_truth/R_reco` against the pre-fix baseline in `MTC_OVERCURVATURE_FINDINGS.md`. If MTC quality is unchanged or improved, the chi² rise is just bookkeeping; if it degrades, the SiPad pad-quantisation bias may need a measurement-σ inflation.

### Issue 4 — what the ACTS output represents

**Status: not yet investigated this session.** Recorded here so it does not get lost.

**Sketch of the expected answer (to be confirmed by reading `ACTSProtoTracker` write-out path and `EDM4HEP2RNTuple` in job5).** The CKF produces a smoothed `Acts::TrackContainer` whose track states carry both the original source links (pointers back to the measurements) and the fitted local parameters at each surface (analogous to the "predicted line" in linear regression — i.e. *new* points, not the raw hits). When `ACTSProtoTracker` writes the EDM4HEP `Track`, it has the choice of emitting the measurement positions, the smoothed-state positions, or both. The next step is to read the write-out loop and document which of these ends up in `ShipHits.root` and under what field name.

---

## See also
- `MTC_OVERCURVATURE_FINDINGS.md` — 2026-05-17 surface-material-attachment investigation (resolved the 3–10× momentum bias).
- `known_issues.md` — running list of remaining quality issues.
- `remaining_acts_issues.md` — the four-issue checklist this section addresses.

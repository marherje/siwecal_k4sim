# ACTS Track Reconstruction Integration

## Overview

This repo reconstructs tracks through the **15 silicon layers of the SiW-ECAL
prototype**, using them as tracking surfaces. There is one detector, `SiPad`,
and no separate tracker: `simulation/geometry/SND_compact.xml` includes only
`SiPadDetector.xml`.

ACTS version: 44.3.0 (from the key4hep stack pinned in `.key4hep-release`).
`CMakeLists.txt` requires `Acts` (Core + PluginDD4hep) and `k4ActsTracking`, so
the whole plugin fails to build without them — including the calorimetry-only
algorithms.

Three components, all in `gaudi_source/`:

| Component | Kind | Role |
|---|---|---|
| `ACTSGeoSvc` | Gaudi service | Builds the `Acts::TrackingGeometry`: one plane surface per ECAL layer, with that layer's material |
| `ShowerTagger` | Gaudi algorithm | Identifies EM cascades and flags their hits so they never become measurements; writes the showers |
| `SiPadMeasConverter` | Gaudi algorithm | Pad hits → `edm4hep::TrackerHit3D` measurements bound to those surfaces, minus the flagged ones |
| `ACTSProtoTracker` | Gaudi algorithm | Hough seeding → CKF → KalmanFitter refit → `edm4hep::TrackCollection` |

Configured by the single shared job `gaudi_jobs/pid2026_common/job4_tracking.py`.
See `docs/gaudi_pipeline.md` for the job's environment variables, property
values and measured performance; this document covers the internals.

> **Heritage note.** This code is a fork of `key4ship_PoC`, which tracked three
> subdetectors (SiTarget strips, SiPad pads, MTC SciFi stereo layers). Anything
> in the sources that mentions stations, stereo angles, planes or U/V pairing is
> vestigial: here `detID` is always 1, `station` and `plane` are always −1, and
> every measurement is a 2D pad hit. The `SiTarget`/`MTC` code paths are gone.

---

## ACTSGeoSvc

**Files:** `gaudi_source/ACTSGeoSvc.{h,cpp}`, interface `gaudi_source/ISNDGeoSvc.h`
(which extends `k4ActsTracking/IActsGeoSvc.h`).

### Configuration

| Property | Description |
|---|---|
| `CompactFile` | Path to the DD4hep compact XML. Pass an absolute path — the service calls `fromXML` directly and the job may run from any directory |

### Interface methods

| Method | Returns |
|---|---|
| `trackingGeometry()` | The `Acts::TrackingGeometry` |
| `allSurfaces()` | The 15 sensitive surfaces, sorted by ACTS X (the beam axis) |
| `surfaceByAddress(detID, station, layer, plane)` | Surface lookup by address; SiPad uses `(1, -1, layer, -1)` |
| `geometryContext()` | The `Acts::GeometryContext` |

### Finding the sensitive surfaces

The service walks the **DD4hep `DetElement` tree**, not TGeo volume names:
`detector("SiPad")` → one child per layer (its `id()` is the layer number) →
one child per slice. A slice is that layer's tracking surface when its volume
**or any descendant** carries the sensitive-detector flag.

Both halves of that condition are load-bearing:

* With the tiled segmentation, `buildWafers` puts the sensitive flag on the
  `<slice>_wafer_pads` volume; the slice container itself is gap material (air).
  A non-recursive check finds nothing.
* The sensitive slice index is **not constant across layers** — 5 for layers
  1–14, 10 for layer 0, which has the extra Al/air entrance slices.

Positions come from `DetElement::nominal().worldTransformation()` and half-sizes
from the slice's `TGeoBBox`; the layer index comes from the DetElement id. No
string parsing anywhere.

If a layer yields anything other than exactly one sensitive slice, the service
fails initialization rather than continuing with a wrong geometry.

> **What this replaced, and why it matters.** The original implementation
> selected surfaces by matching volume names against a hardcoded `"_slice_4"`.
> When the segmentation changed (July 2026) the silicon moved to slice 5/10 and
> `_slice_4` became an *air* slice — so the service still found exactly 15
> surfaces (the count matched, which is why nobody noticed) but every one of
> them sat at the wrong z. Nothing errored. Asking DD4hep for sensitivity makes
> that class of failure impossible.

### Geometry construction

The builder is manual rather than `Acts::convertDD4hepDetector` (see
`docs/acts_material_migration.md`, Path A, still an open follow-up). Pieces that
look redundant but are not:

| Piece | Why |
|---|---|
| `SNDDetectorElement` | The CKF classifies any surface with `associatedDetectorElement() == nullptr` as **passive**. Constructing the `PlaneSurface` with a detector element is what makes the surfaces sensitive |
| `SurfaceArray` on every `PlaneLayer` | Without a non-null surface array, `Layer::resolve()` always returns false and the CKF visits zero surfaces |
| `rot90Y = R_Y(pi/2)` | Maps the surface normal onto the beam so tracks run along ACTS X, avoiding the theta = 0 singularity of ACTS bound coordinates. It also maps local x to **minus** global Z — see the sign warning below |
| Two `Acts::navigation` layers | At the z extremes, for correct extrapolation |
| A single `CuboidVolumeBounds` `TrackingVolume` | Sidesteps a `CuboidVolumeBuilder` position-offset bug |
| `LayerArrayCreator` with `Acts::arbitrary` binning on `AxisX` | Respects the actual (non-uniform) layer positions |

### Surface material

Each surface carries its layer's **entire slice stack** — tungsten absorber,
silicon, PCB, Cu, carbon fibre; air and vacuum skipped — folded into one
`Acts::HomogeneousSurfaceMaterial` via `MaterialSlab::combineLayers`.
Thicknesses come from the slice `TGeoBBox` shapes and X0/L0/A/Z/rho from
`dd4hep::Detector::material()`, so nothing about the material is hardcoded and
the budget tracks the XML.

| Layers | t | X0 | t/X0 |
|---|---|---|---|
| 0 | 14.4 mm | 7.06 mm | 2.04 |
| 1–7 | 8.1 mm | 6.64 mm | 1.22 |
| 8–13 | 9.5 mm | 5.88 mm | 1.62 |
| 14 | 11.65 mm | 7.10 mm | 1.64 |

This matters because the CKF's multiple-scattering and energy-loss flags read
`surface.surfaceMaterial()`; with no material attached they are silent no-ops.
The service logs the table above at INFO on every run — check it after any
geometry change.

### Unit notes

DD4hep and TGeo work in cm, ACTS in mm: every length is multiplied by 10 on the
way in. `radLength()`/`intLength()` are likewise cm.

---

## Coordinates — the one thing to get right

Two frames are in play, and the mapping between them is **not** a plain axis
relabelling.

* **DD4hep / detector:** z is the beam, (x, y) transverse.
* **ACTS global:** X is the beam (after `rot90Y`), so the transverse plane is
  (Y, Z).

`rot90Y = R_Y(pi/2)` sends the local axes to

```
local x -> (0, 0, -1)     local y -> (0, 1, 0)     local z -> (1, 0, 0)
```

so local z is the beam as intended, but **local x maps to minus global Z**.
Inverting: a global offset `(0, dy, dz)` from the surface centre has bound
coordinates `(loc0, loc1) = (-dz, dy)`.

`SNDCalibrator` defines the measurement as `(loc0, loc1) = (dd_x, dd_y)`.
Consistency therefore requires

```
global Y = +dd_y        global Z = -dd_x
```

for both the seed position and the seed direction. Filling `ePos2` with `+dd_x`
puts the seed `2*|dd_x|` away from its own track — invisible on the beam axis,
but 88 mm for a muon at `dd_x = -44` mm, whose first-layer chi2 `(88/10.1)^2 = 76`
then exceeds `Chi2CutOff = 70`, so the CKF loses the track completely. This cost
~9% of the muon efficiency and was invisible in any single-event display near
the beam centre.

On output the convention is inverted back, so `TrackState::referencePoint`
carries DD4hep `(x, y, z_beam)`.

---

## ShowerTagger — why an EM event yields no track

**File:** `gaudi_source/ShowerTagger.cpp`

A track through an electromagnetic cascade is not a physical object. Once a
particle showers, the pads it lights are secondaries spraying transversely, not
samples of a trajectory; a line fitted through them is a well-formed track with
a small chi2 and no meaning. The physical answer for an EM event is "this is a
shower, with this energy and barycentre".

In a SiW-ECAL that is not an edge case. Every layer is 1.2-2.0 X0, so a
high-energy electron is already showering in layer 0 or 1 and there is simply no
incoming segment to fit. Filtering shower hits *inside* the tracker is too late
and too weak — density cuts trim the core but the remaining fringe still
supports several plausible "tracks". The hits have to be kept out of the sampling
altogether.

**Onset.** Hits are counted per layer. A MIP lights one or two pads; a shower
lights tens. The onset is the first run of `ShowerMinConsecutive` (2) consecutive
layers with at least `ShowerNHitsThreshold` (4) hits — consecutive in layer
*number*, so a gap breaks the run, and two layers so a single delta-ray spike
does not trigger it.

**Veto.** Hits from the onset onwards are flagged. If fewer than
`MinTrackLayers` (4) layers precede the onset, every hit is flagged: one to three
points cannot define a trajectory, and passing them on is exactly how a shower
event acquires a "track".

| Particle | Outcome |
|---|---|
| Muon | No onset; the fifteen-layer track is untouched |
| Electron | Onset at layer 0-1 → one shower, **zero tracks** |
| Radiative muon | Late onset → incoming track kept *and* a shower |
| Pion punching through then showering | Same: incoming track plus shower — what PID wants |

**Outputs.** `OutputFlags` is a `podio::UserDataCollection<int32_t>`, one entry
per input hit, 0 = track-like and 1 = shower — the same idiom as
`ChannelMapper`'s `OutputMaskedFlags`. It reaches ACTS through
`SiPadMeasConverter.InputFlags`; because the flags are positional, the converter
fails loudly on a length mismatch rather than vetoing the wrong hits.

`OutputShowers` is an `edm4hep::ClusterCollection`: `type = 1`, `energy` = total
in the input's units (MIPs after digitisation), `position` = energy-weighted
barycentre, and `shapeParameters` = [start layer, layer of maximum, layers
spanned, transverse RMS in mm, number of hits]. Hit relations are not filled:
`Cluster::addToHits` wants `CalorimeterHit` and the input here is
`SimCalorimeterHit`.

`Enabled = False` turns tagging off entirely (all flags zero, no showers), which
is the way to reproduce the old behaviour for a comparison.

Measured on a 300-event e- 74 GeV slice: 299 showers, **0 tracks**, median 4990
MIP with the onset at layer 1 and the maximum at layer 6. On 1000 mu- events:
49 showers (radiative muons), 987 tracks — the 13 losses are muons that dumped a
cascade in the first few layers.

The tracker's own `HitPurgeWindow` / `IsolationWindow` filters remain as
fine-grained cleanup for delta rays inside an otherwise track-like event; with
`ShowerTagger` upstream they no longer carry the shower problem on their own.

---

## SiPadMeasConverter

**File:** `gaudi_source/SiPadMeasConverter.cpp`

Reads a pad-hit collection and writes `edm4hep::TrackerHit3D`:

* `layer` is decoded from the cell ID with a `dd4hep::DDSegmentation::BitFieldCoder`
  built from the `BitField` property, and written into `quality` — that is how
  `ACTSProtoTracker` finds the surface, via `surfaceByAddress(1, -1, layer, -1)`.
* Position is taken **from the hit** (`getPosition()`), not recomputed from the
  cell ID.
* Covariance is `pitch²/12` on each transverse axis, from `PixelSizeX/Y`.

`InputFlags` names the per-hit veto collection (`ShowerTagger`'s
`OutputFlags`); flagged hits are counted as `vetoed(shower)` in the DEBUG line
and never become measurements. Leaving it empty disables the veto.

`BitField` and `PixelSize*` come from `parse_geometry` in the job, i.e. straight
out of the compact XML's `<readout>` block, so they cannot drift from the
segmentation the way the previously hardcoded copies did.

**Input collection:** must be a *pre-flip* collection (`SiPadHitsWindowed` or
`SiPadHitsDigi`). `DetectorFlipper` rewrites the hit z into the test-beam frame,
which no longer matches the ACTS surfaces.

---

## ACTSProtoTracker

**File:** `gaudi_source/ACTSProtoTracker.cpp`

### Internal types

| Type | Role |
|---|---|
| `SNDMeasurement` | A 2D pad measurement: surface, `(localCoord, localCoord2)`, variances, time, energy |
| `SNDSourceLink` | Index into the measurement vector plus the surface geoID |
| `SNDSourceLinkAccessor` | `lower_bound`/`upper_bound` over the geoID-sorted source links: surface → measurement range in O(log N) |
| `SNDFixedNavigator` | Wraps `Acts::DirectNavigator`, injecting the surface list at `makeState()` |
| `SNDCalibrator` | Writes the calibrated 2D coordinates, covariance and projector subspace |
| `SNDSurfaceAccessor` | Source link → surface; required by `KalmanFitterExtensions` |
| `IronSlabBField` | `MagneticFieldProvider` kept from the SND setup; with `IronFieldRanges` empty the field is zero, which is correct for a test beam |

### Two non-obvious obligations

**`SNDFixedNavigator::endOfWorldReached()` must return `navigationBreak`.**
`Propagator::propagate()` leaves its stepping loop only when an aborter fires,
and the CKF's only geometric aborter is `EndOfWorldReached`, which calls exactly
this method. For a fixed surface sequence, running out of surfaces *is* the end
of the world. Returning a hardcoded `false` leaves the propagator free-stepping
at `maxStepSize` after `DirectNavigator` has set `navigationBreak`, until
`maxSteps` is reached and `PropagatorError::StepCountLimitReached` comes back —
the CKF then discards the branch and **every event produces zero tracks**. (In
the PoC, with 150 surfaces over 4 m, this surfaced as an intermittent
`PropagatorError:2` on ~6% of events and was never closed; with only 15 surfaces
it is deterministic.)

**`SNDCalibrator` must call `ts.setUncalibratedSourceLink()`.**
`Acts::TrackStateCreator` leaves that to the calibrator by convention. Without
it, `hasUncalibratedSourceLink()` is always false, so every hit fingerprint is
empty — which silently disables duplicate rejection (`overlapFraction` is always
0) and the frozen-hit refit.

### Algorithm flow

1. **Collect measurements** from the input collection; surfaces resolved by
   address from the layer index. Sorted by surface position along the beam.
2. **Shower-hit purge** (`HitPurgeWindow` / `HitPurgeMaxNeighbors`): per surface,
   drop a hit when more than `HitPurgeMaxNeighbors` other hits lie within the
   window. Unlike the seed-level isolation filter this removes hits from the
   *fit*, so a shower core cannot contaminate a candidate.
3. **Hough seeding** (`findSeeds`): histogram the pad hits in (x, y), take local
   maxima with non-maximum suppression, reject peaks whose crossings-per-layer
   multiplicity marks them as showers, then fit a straight line through each
   peak's compatible hits to get the entry point on the first surface **and** a
   direction. (The direction used to be hardcoded along the beam, which is right
   only for a track exactly parallel to it.)
4. **CKF** (`runCKFPass`) per seed: `EigenStepper` + `SNDFixedNavigator`,
   `GainMatrixUpdater`, `MeasurementSelector(Chi2CutOff, NumMeasCutOff)`,
   multiple scattering and energy loss on. Keeps the candidate with the **best**
   chi2/ndf that has at least 3 measurements.
5. **Final refit** (`FinalRefit`): refit the frozen hit set with
   `Acts::KalmanFitter` + `GainMatrixSmoother`, re-seeded from the smoothed CKF
   state, for unbiased parameters and covariances. Falls back to the CKF track.
6. **Acceptance** on `MaxChi2PerNdf`, with `ndf = sum(calibratedSize) - 5`.
7. **Seed cleaning** (`SeedCleaning`): drop the accepted track's hits from the
   source-link pool so later seeds cannot re-find the same particle.
8. **Deduplication** at end of event, best-first: candidates sorted by chi2/ndf
   ascending, a candidate dropped when it shares more than
   `DuplicateOverlapFraction` of the smaller hit set with an already-kept one.
   The survivor of each group is therefore always the best fit, independent of
   seed order.

Outlier states need explicit handling when building a hit fingerprint:
`TrackStateCreator` sets `MeasurementFlag` on outliers **as well as**
`OutlierFlag`, and outliers must not reach the frozen-hit refit.

### Neighbour windows are in pad pitches

`IsolationWindow` and `HitPurgeWindow` count neighbours on a plane, so they must
be expressed in units of the pad pitch. **Any window ≤ one pitch (5.53 mm)
counts zero neighbours by construction** — the nearest other pad is exactly one
pitch away — which silently disables the filter. The 5.0 mm value inherited from
the PoC came from the SND SiTarget, where the strips are 75 µm and 5 mm spans
dozens of channels. The job uses 1.5 pitches, reaching the 4 orthogonal (1 pitch)
and 4 diagonal (1.41 pitches) neighbours, i.e. the 8 pads around a hit.

### Track output

One `edm4hep::Track` per surviving candidate, with `type = 1`, `chi2`, and `ndf`
already reduced by the 5 helix parameters (so consumers divide directly):

| TrackState | Contents |
|---|---|
| `AtIP` (one) | `D0`, `Z0` = the seed's transverse `(x, y)` |
| `AtOther` (one per surface) | `referencePoint` = DD4hep `(x, y, z_beam)`; `phi`, `tanLambda` from the smoothed state; `D0` = raw `loc0`; `Z0` = that state's innovation chi2 |

`Z0` on `AtOther` carrying a per-state chi2 is deliberate: the track total is
useless for finding *which* surface drives the fit quality.

### Debug output

`TRACKING_LOGLEVEL=DEBUG` enables the per-seed CKF summary, the Hough peak list
with multiplicities, the `DIAG-KF` refit line, and `DIAG-TS`, which dumps every
track state of the first three events — measurement/outlier/hole flags,
predicted vs measured local coordinates and predicted sigmas. That dump is the
fastest way to tell a rejected-by-chi2 hit from a genuinely missing one.

---

## Interpreting chi2/ndf on this detector

The chi2/ndf median sits near **zero**, and that is geometry, not a bug: with a
5.53 mm pitch and 15 mm between layers, a track must be tilted by more than
`atan(5.53/15)` ≈ 20° to change pad. A beam muon therefore deposits in the *same*
pad on all 15 layers, the measurements are literally identical, and a straight
line passes exactly through them. Only genuinely tilted tracks populate the tail.

Consequence: `MaxChi2PerNdf` is a weak discriminant here, and the useful quality
handles are the number of measurements, the number of holes, and the residual
against the hits.

---

## See also

* `docs/gaudi_pipeline.md` — the job, its properties and measured performance
* `docs/acts_material_migration.md` — Path B (implemented) vs Path A
  (`Acts::convertDD4hepDetector`, still open)
* `docs/dd4hep_plugins.md` — the `SiPadDetector` plugin and its wafer tiling

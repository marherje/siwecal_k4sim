# Gaudi Reconstruction Pipeline

## Overview

Sequential `k4run` jobs turn ddsim output into digitised hits, tracks and
analysis trees. Only `SiPad` (the SiW-ECAL prototype) exists in this repo.

```
output_*.edm4hep.root
  → job1: EventShuffler       → shuffled.edm4hep.root        (optional: pile-up / multi-source)
  → job2: EventWindowSplitter → timewindows.edm4hep.root     (optional: time windows)
  → job3: GeV2MIPConversion + BasicDigitizer + DetectorFlipper + ChannelMapper
                              → digitized.edm4hep.root
  → job4: SiPadMeasConverter + ACTSProtoTracker → tracks.edm4hep.root
  → job5: EDM4HEP2RNTuple     → ShipHits.root                (optional)
```

Not every pipeline runs every job. The beam and per-chunk productions go
straight from ddsim to job3, then job4, then the ecal-tree conversion; jobs 1, 2
and 5 are used by the `1_*_PG*` pipelines. Jobs 3 and 4 are single shared
configs under `gaudi_jobs/pid2026_common/`; see the per-pipeline `.sh` scripts
for the exact chain.

```bash
# muon beam, end to end
bash gaudi_jobs/1_mu_beam_pipeline/1_mu_pipeline.sh
```

---

## Job 1 — EventShuffler

**File:** `gaudi_source/EventShuffler.cpp`  
**Config:** `gaudi_jobs/*/job1_shuffler.py`

Merges N simulation files into one super-event. Assigns source IDs and time offsets per file to simulate pile-up.

### Architecture note
All work happens in `finalize()`. The algorithm reads files directly via `podio::ROOTReader`, bypassing Gaudi's IOSvc entirely. `execute()` is a no-op.

```python
ApplicationMgr(TopAlg=[shuffler], ExtSvc=[], EvtSel="NONE", EvtMax=1)
# No IOSvc — would conflict with direct podio I/O
```

### Key properties
| Property | Description |
|----------|-------------|
| `InputFiles` | List of edm4hep ROOT files (one per source) |
| `SourceIDs` | Integer ID for each file (same order) |
| `Delays` | Inter-event time delay in ns for each source |
| `CollectionsSiPad` | Collection name for each input file |
| `OutputFile` | Output file path |
| `MaxEventsPerSource` | 0 = no limit |

### Source ID encoding
`edm4hep::CaloHitContribution` has no source field. Source ID is stored in the **PDG field** (`contrib.setPDG(source_id)`). PDG is unused elsewhere in this pipeline.

---

## Job 2 — EventWindowSplitter

**File:** `gaudi_source/EventWindowSplitter.cpp`  
**Config:** `gaudi_jobs/*/job2_splitter.py`

Splits the merged super-event into 25 ns time windows. Each window becomes one EDM4HEP frame in the output file. Uses Gaudi IOSvc (normal I/O, not bypass mode).

### Key properties
| Property | Description |
|----------|-------------|
| `WindowSize` | Time window width in ns (default 25) |
| Input/output collections | Configured via IOSvc `keep` rules |

---

## Job 3 — Digitization

**Files:** `gaudi_source/GeV2MIPConversion.cpp`, `BasicDigitizer.cpp`,
`DetectorFlipper.cpp`, `ChannelMapper.cpp`
**Config:** `gaudi_jobs/pid2026_common/job3_digitize.py` (shared)

```
SiPadHits → GeV2MIPConversion → SiPadHitsMIP
          → BasicDigitizer    → SiPadHitsDigi      <-- tracking input (pre-flip)
          → DetectorFlipper   → SiPadHitsFlipped
          → ChannelMapper     → SiPadHitsMapped + SiPadHitsMasked
```

| Algorithm | Role |
|---|---|
| `GeV2MIPConversion` | Energy in GeV → MIPs. `MIPValues` takes one value per layer (from `mip_extraction_pipeline`); `MIPValue` is the scalar fallback |
| `BasicDigitizer` | Applies the MIP `Threshold`, dropping hits below it |
| `DetectorFlipper` | Rewrites hit z into the **test-beam frame** using `mappings/slab_z_positions.yml` — the single source of truth for the per-slab z, shared with the event viewer |
| `ChannelMapper` | Cell IDs → test-beam format (`system:8,slab:8,chip:16,channel:8,sca:8`) via the pad maps in `mappings/`, and masks dead channels from the muon calibration tree (`CalibThreshold`, e.g. `th230`, masks ~3.5% of channels) |

**Tracking reads `SiPadHitsDigi`, before `DetectorFlipper`**: the flip moves the
hit z into the test-beam frame, which no longer matches the ACTS surfaces built
from the compact XML.

---

## Job 4 — Tracking

**Files:** `gaudi_source/SiPadMeasConverter.cpp`, `gaudi_source/ACTSProtoTracker.cpp`, `gaudi_source/ACTSGeoSvc.cpp`
**Config:** `gaudi_jobs/pid2026_common/job4_tracking.py` — a single job shared by every pipeline

There is exactly one tracking job. It used to be copy-pasted into each
`gaudi_jobs/1_*_pipeline/` directory with the bit field and pad pitch hardcoded;
those copies went stale as soon as the segmentation changed. Everything that
describes the detector now comes from `simulation/geometry/parse_geometry.py`,
and everything that varies per pipeline comes from environment variables:

| Variable | Default | Meaning |
|---|---|---|
| `INPUT_FILE` | `timewindows.edm4hep.root` | input edm4hep file |
| `OUTPUT_FILE` | `tracks.edm4hep.root` | output edm4hep file |
| `INPUT_COLLECTION` | `SiPadHitsWindowed` | hit collection to track on |
| `SEED_MOMENTUM` | `3.0` | beam momentum [GeV] |
| `TRACKING_LOGLEVEL` | `INFO` | `DEBUG` for the per-surface dumps |

**Which input collection:** pipelines that split time windows track on
`SiPadHitsWindowed`; the per-chunk production tracks on **`SiPadHitsDigi`**.
Always a *pre-flip* collection: `DetectorFlipper` rewrites the hit z into the
test-beam frame, which no longer matches the ACTS surfaces — those come from
the same compact XML as the simulation.

**Step 1 — `SiPadMeasConverter`:** turns the pad hits into
`edm4hep::TrackerHit3D` measurements, writing the layer index into `quality` so
the tracker can look the surface up by address. Variance is `pitch²/12` per
axis. Positions are taken from the hit, not recomputed.

**Step 2 — `ACTSProtoTracker`:** Hough seeding → CKF → KalmanFitter refit →
event-level deduplication, writing `ACTSTracks` with one `AtIP` state carrying
the seed position and one `AtOther` state per surface.

### Geometry: one surface per layer, from DD4hep

`ACTSGeoSvc` walks the **DD4hep `DetElement` tree** (`SiPad` → layer → slice)
and picks the slice whose volume *or any descendant* is sensitive. Both halves
matter: the sensitive flag sits on the `*_wafer_pads` volume created by
`buildWafers`, not on the slice container, and the sensitive slice index is not
even constant across layers (5 for layers 1–14, 10 for layer 0). The previous
implementation matched volume **names** against a hardcoded `_slice_4`; when the
segmentation changed it silently returned 15 surfaces sitting on *air* slices at
the wrong z, with no error anywhere.

Each surface carries that layer's **whole slice stack** (W, Si, PCB, Cu, CF —
air skipped) combined into one `HomogeneousSurfaceMaterial`, with thicknesses
read from the `TGeoBBox` shapes and X0/L0/A/Z/rho from DD4hep. Result:
t/X0 = 1.22 for layers 1–7, 1.62 for 8–13, 2.04 for layer 0.

### CKF architecture

| Component | Location | Role |
|-----------|----------|------|
| `SNDFixedNavigator` | `ACTSProtoTracker.cpp` | Wraps `DirectNavigator`; injects the 15-surface list at `makeState()` so the CKF's `setPlainOptions()` cannot erase it. **`endOfWorldReached()` must report `navigationBreak`** — see below |
| `SNDDetectorElement` | `ACTSGeoSvc.cpp` | `DetectorElementBase` subclass; makes `associatedDetectorElement() != nullptr` so `CKFActor` treats surfaces as sensitive instead of passive |
| `SNDSourceLinkAccessor` | `ACTSProtoTracker.cpp` | Binary-search lookup: surface geoID → measurement range in O(log N) |
| `SNDCalibrator` | `ACTSProtoTracker.cpp` | Sets the calibrated 2D coordinates and the projector subspace. **Must call `setUncalibratedSourceLink()`** — `TrackStateCreator` leaves that to the calibrator, and without it every hit fingerprint is empty, which silently disables duplicate rejection and the frozen-hit refit |
| `SNDSurfaceAccessor` | `ACTSProtoTracker.cpp` | Resolves a source link back to its surface; required by `KalmanFitterExtensions` for the final refit |
| `IronSlabBField` | `ACTSProtoTracker.cpp` | Kept from the SND setup; with `IronFieldRanges` empty (the test-beam case) the field is a zero `ConstantBField` |

Three things are easy to get wrong here and produce **no error at all**:

1. **`endOfWorldReached()`.** `Propagator::propagate()` only leaves its stepping
   loop when an aborter fires, and the CKF's only geometric aborter is
   `EndOfWorldReached`, which calls `navigator.endOfWorldReached()`. For a fixed
   surface sequence, running out of surfaces *is* the end of the world. Returning
   a hardcoded `false` leaves the propagator free-stepping at `maxStepSize` after
   `DirectNavigator` has set `navigationBreak`, until it hits `maxSteps` and
   returns `PropagatorError::StepCountLimitReached` — every event yields zero
   tracks.
2. **The seed coordinate mapping.** The geometry applies `rot90Y = R_Y(pi/2)`,
   which maps the local x axis to **minus** global Z. Since the calibrator
   defines the measurement as `(loc0, loc1) = (dd_x, dd_y)`, consistency
   requires `global Z = -dd_x` and `global Y = dd_y`. Getting the sign wrong puts
   the seed `2*|dd_x|` away from its own track: harmless on the beam axis, but
   88 mm for a muon at `dd_x = -44` mm, whose first-layer chi2 then exceeds
   `Chi2CutOff` and the track is lost. The same mapping applies to the direction.
3. **Neighbour-window units.** `IsolationWindow` and `HitPurgeWindow` count
   neighbours on a plane, so they must be expressed in **pad pitches**. Any
   window <= one pitch (5.53 mm) counts zero neighbours by construction and
   silently disables the filter. The job uses 1.5 pitches, which reaches the 8
   pads surrounding a hit.

### ACTSProtoTracker key properties

| Property | Job value | Description |
|----------|-----------|-------------|
| `AutoSeed` | `True` | Hough-transform seeding |
| `MaxSeeds` | 3 | Maximum seeds tried per event |
| `HoughBinSize` | 5.0 | Hough accumulator bin size [mm] |
| `HoughHalfSize` | `geo.hough_half_size` (95.53) | Accumulator range [mm], from the transverse envelope |
| `HoughMinVotes` | 3 | Minimum votes to form a seed |
| `SeedCompatRadius` | 8.0 | Radius [mm] for hits compatible with a peak |
| `SeedStripPitch` | pad pitch | Bin for the seed-position refinement |
| `SeedMomentum` | `$SEED_MOMENTUM` | Seed momentum [GeV] |
| `Chi2CutOff` | 70.0 | `MeasurementSelector` per-surface chi2 cut |
| `NumMeasCutOff` | 1 | Max measurements accepted per surface |
| `MaxChi2PerNdf` | 10.0 | Track acceptance on chi2/ndf, ndf = sum(calibratedSize) - 5 |
| `HoughMaxMultiplicity` | 10.0 | Max crossings/layer for a peak to be a track candidate |
| `IsolationWindow` / `IsolationMaxNeighbors` | 1.5 pitch / 2 | Seed-level shower rejection |
| `HitPurgeWindow` / `HitPurgeMaxNeighbors` | 1.5 pitch / 4 | Pool-level shower rejection: drops a hit when more than half of the 8 pads around it are lit |
| `SeedCleaning` | `True` | Removes an accepted track's hits from the pool before the next seed |
| `FinalRefit` | `True` | Refits the frozen hit set with `Acts::KalmanFitter` for unbiased parameters |
| `DuplicateOverlapFraction` | 0.7 | Event-level best-first deduplication threshold |
| `IronFieldRanges` | `[]` | Per-slab field map; empty = no field, correct for a test beam |

### Seeding

The Hough transform histograms the 2D pad hits in (x, y) and takes local
maxima. For each peak, a straight line is fitted through the compatible hits to
give the seed both its entry point on the first surface **and** a direction —
the seed direction used to be hardcoded along the beam, which only works for a
track exactly parallel to it.

### Performance

Measured on `1_mu_beam_pipeline` (1000 events, mu- 100 GeV) and a 300-event
slice of an e- 74 GeV chunk:

| | mu- 100 GeV | e- 74 GeV |
|---|---|---|
| Events with a track | 999/1000 | 298/300 |
| Tracks per event | 1.00 | 2.11 |
| Track spans all 15 layers | 99.8% | — |
| Residual vs hits (rms) | 0.21 mm (x), 0.28 mm (y) | — |
| Best track within one pad of the layer-0 hit | — | 95.6% |
| Hits purged as shower-like | ~0% | 30.5% |

chi2/ndf has a median near zero and that is expected, not a bug: with a 5.53 mm
pitch and 15 mm between layers a track must be tilted by more than
atan(5.53/15) ~ 20 deg to change pad, so a beam muon gives *identical*
measurements on all 15 layers and the residual is zero by construction. The
tail comes from genuinely tilted tracks.

### ACTSGeoSvc in job config
```python
from Configurables import ACTSGeoSvc
geo = ACTSGeoSvc("ACTSGeoSvc")
geo.CompactFile = str(COMPACT_FILE)      # absolute: the job runs from anywhere
ApplicationMgr(..., ExtSvc=[iosvc, geo])  # Service in ExtSvc, not TopAlg
```

---

## Job 5 — EDM4HEP2RNTuple

**File:** `gaudi_source/EDM4HEP2RNTuple.cpp`  
**Config:** `gaudi_jobs/*/job5_rntuple.py`

Converts EDM4HEP collections to a ROOT RNTuple (`ShipHits.root`) for analysis.

### Written collections
- `SiPad` — from `SiPadHitsWindowed` (`Collections` / `BitFields`)
- measurements — from `SiPadMeasurements` (`MeasCollections` / `MeasBitFields`)
- `Tracks` — from `ACTSTracks` in `TrackFile` (`TrackCollectionName`)

Only the `1_*_PG*` pipelines run this job; the beam and per-chunk productions
use `analysis/sim_to_ecal_tree` on `SiPadHitsMapped` instead.

---

## Adding a New Gaudi Algorithm

1. Create `gaudi_source/MyAlgorithm.cpp` with `DECLARE_COMPONENT(MyAlgorithm)` at end
2. Add `.cpp` to `gaudi_add_module(SND_reco SOURCES ...)` in `CMakeLists.txt`
3. `/build` to rebuild
4. Import in Python: `from Configurables import MyAlgorithm`

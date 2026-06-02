# ACTS material-budget migration plan

**Status:** Path B implemented 2026-05-17. Path A remains a future follow-up.
Drafted and implemented in response to the over-curvature investigation
([`gaudi_jobs/muon_pipeline/MTC_OVERCURVATURE_FINDINGS.md`](../gaudi_jobs/muon_pipeline/MTC_OVERCURVATURE_FINDINGS.md))
and the latent-risk note in
[`docs/known_issues.md`](known_issues.md).

**Path B result:** median |R_truth/R_reco| = 0.94 (was 3–10×); median reco |dx|
through MTC = 1–6 mm (was 1600–3400 mm). See `known_issues.md` for full table.

## TL;DR

Right now `ACTSGeoSvc` builds the ACTS `TrackingGeometry` from a manual TGeo
walk and attaches **no material** to any surface or layer. As a consequence,
the CKF's MS/EL flags
(`gaudi_source/ACTSProtoTracker.cpp:1340-1341`,
`/*multipleScattering=*/true, /*energyLoss=*/true`) are no-ops — the fitter
treats 50 mm of MTC iron per layer as vacuum, absorbs each MS kick into the
fitted q/p, and reconstructs a track with ~3-10× too tight a radius. A blanket
"iron everywhere" patch is wrong (SiTarget is silicon, SiPad has tungsten +
silicon, only MTC has iron). The right fix is to make the ACTS geometry
materially consistent with what Geant4/DD4hep already simulates.

Two paths:

- **Path B (recommended near-term):** keep the current geometry builder; query
  DD4hep for each plane's parent volume material and attach
  `HomogeneousSurfaceMaterial` directly. Scoped, reversible, doesn't disturb
  the workarounds for ACTS' `CuboidVolumeBuilder` offset bug or the
  `SNDFixedNavigator`.

- **Path A (recommended follow-up):** annotate the three detector plugins with
  `Acts::ActsExtension` at the right `DetElement` granularity and replace
  `ACTSGeoSvc::initialize` with a single call to
  `Acts::convertDD4hepDetector(...)`. One source of truth for both simulation
  and reconstruction geometry. Larger blast radius, larger payoff.

This document specifies both paths in enough detail to implement, with
verification checkpoints at each step.

---

## Background — current state

`gaudi_source/ACTSGeoSvc.cpp` (current):

1. Loads DD4hep geometry from the compact file
   (`gaudi_source/ACTSGeoSvc.cpp:80-82`).
2. Walks the TGeo tree, extracts sensitive plane positions and half-sizes into
   a `PlaneInfo` list (`gaudi_source/ACTSGeoSvc.cpp:86-238`).
3. Builds one `Acts::PlaneLayer` per plane via
   `SNDDetectorElement` + `PlaneSurface(bounds, detElement)`
   (`gaudi_source/ACTSGeoSvc.cpp:263-303`).
4. Builds a single `Acts::CuboidVolumeBounds` `TrackingVolume` around all
   layers (`gaudi_source/ACTSGeoSvc.cpp:341-369`).
5. Returns this geometry via `ISNDGeoSvc::trackingGeometry()` and the surface
   list via `allSurfaces()`/`surfaceByAddress()`.

What's missing:

- No `Acts::ISurfaceMaterial` is ever attached to any surface or layer.
- No passive layers for the inter-plane material (50 mm Fe between MTC SciFi
  planes; W absorber between SiPad planes).
- `Acts::TrackingVolume` has no `IVolumeMaterial`.

CKF MS noise is computed in
`Acts/include/Acts/TrackFitting/detail/GsfActor.hpp` (and equivalent in CKF)
by reading `surface.surfaceMaterial()`. If null, the contribution is zero.

---

## Material budgets per subdetector

Pulled from `simulation/geometry/*.xml` and the compact include hierarchy
(`SND_compact.xml`). Numbers below are the design intent; the exact thicknesses
should be confirmed in the XML before encoding any constants.

| Subdetector | Compact XML | Sensitive plane material | Sensitive thickness | Upstream material (per layer) |
|---|---|---|---|---|
| SiTarget | `SiTarget.xml` | Silicon | ~0.3 mm | none (just air) |
| SiPad | `SiPadDetector.xml` | Silicon | ~0.5 mm | Tungsten absorber, ~3.5 mm |
| MTC SciFi | `MTCDetector.xml` | Scintillating fibre (≈polystyrene) | ~0.5 mm | Iron, 50 mm |
| MTC scintillator | `MTCDetector.xml` | Polystyrene scint | ~5 mm | (counted under SciFi) |

The MS-dominant material in each case is the *upstream* layer — for SiTarget
that's nothing, for SiPad it's the tungsten, for MTC it's the iron. The
sensitive layer itself is in all three cases thin (≤ 0.5 mm Si or
scintillator) and contributes a small amount of process noise relative to its
upstream absorber.

For ACTS' `Acts::Material` / `Acts::MaterialSlab` constructors we'll need:
- atomic number `Z`
- mass number `A`
- mass density ρ [g/cm³]
- radiation length X₀ [mm]
- nuclear interaction length λ₀ [mm]

DD4hep already has all of these encoded in its `Material` objects. The
recommended source-of-truth lookup is `dd4hep::Material(name).radLength()`
etc. — never hard-code these.

---

## Path B — minimal-touch / material-only graft

### Scope

Modify `gaudi_source/ACTSGeoSvc.cpp` only. No changes to detector plugins,
compact XML, navigator, or volume builder. Add material on each surface after
it's been constructed (after the `auto layer = Acts::PlaneLayer::create(...)`
block at line 294).

### Step 1 — material accessor helper

In `ACTSGeoSvc.cpp`, add a private helper:

```cpp
Acts::MaterialSlab makeSlab(const std::string& dd4hepMaterialName,
                            double thicknessMm) const {
  auto& desc = dd4hep::Detector::getInstance();
  const dd4hep::Material mat = desc.material(dd4hepMaterialName);
  // dd4hep returns lengths in cm; ACTS wants mm.
  const double X0  = mat.radLength()           * 10.0;  // cm → mm
  const double L0  = mat.intLength()           * 10.0;  // cm → mm
  const double rho = mat.density();                     // g/cm³ — ACTS unit OK
  const double A   = mat.A();
  const double Z   = mat.Z();
  return Acts::MaterialSlab(
      Acts::Material::fromMassDensity(X0, L0, A, Z, rho),
      thicknessMm * Acts::UnitConstants::mm);
}
```

(Confirm the unit convention of `dd4hep::Material::radLength()` against the
2026-02-01 release — the conversion may already return cm; the rule is "never
trust units, verify by printing one slab's X₀ in mm and comparing to PDG").

### Step 2 — attach material per surface

After the `allLayers.push_back(layer);` line at `ACTSGeoSvc.cpp:302`, but
*before* committing the loop, attach material based on `pi.detID`:

```cpp
std::shared_ptr<const Acts::ISurfaceMaterial> surfMat;

switch (pi.detID) {
  case 0:  // SiTarget: just the silicon sensor itself
    surfMat = std::make_shared<Acts::HomogeneousSurfaceMaterial>(
        makeSlab("Silicon", /*thickness mm=*/0.3));
    break;

  case 1:  // SiPad: tungsten absorber upstream + silicon sensor
    // Combine into one effective MaterialSlab via thickness-weighted average,
    // OR attach two separate slabs to two surfaces (preferred but invasive).
    // Near-term shortcut: one combined slab on the Si sensor surface that
    // approximates the integrated absorber + sensor pre-step material.
    surfMat = std::make_shared<Acts::HomogeneousSurfaceMaterial>(
        Acts::MaterialSlab::averageLayers(
            makeSlab("Tungsten", 3.5),
            makeSlab("Silicon", 0.5)));
    break;

  case 2:  // MTC SciFi: 50 mm Fe upstream + 0.5 mm scifi
    if (pi.layerInDet == 0) {
      // first SciFi plane in a station: no upstream iron
      surfMat = std::make_shared<Acts::HomogeneousSurfaceMaterial>(
          makeSlab("Polystyrene", 0.5));
    } else {
      surfMat = std::make_shared<Acts::HomogeneousSurfaceMaterial>(
          Acts::MaterialSlab::averageLayers(
              makeSlab("Iron", 50.0),
              makeSlab("Polystyrene", 0.5)));
    }
    break;
}

if (surfMat) {
  detElem->surface().assignSurfaceMaterial(surfMat);
}
```

Notes:
- `assignSurfaceMaterial` is the ACTS API for attaching after construction.
- `MaterialSlab::averageLayers` (or its equivalent in 44.x —
  may be `combine`/`fromLayers`; check the header
  `Acts/Material/MaterialSlab.hpp`) folds two slabs into a single effective
  one. This is the standard ACTS pattern when one wants the pre-step material
  represented on the downstream surface.
- For MTC, attaching the iron to *every* SciFi U and V plane double-counts the
  scattering (each layer has 2 planes 2.35 mm apart but only one 50 mm iron
  slab upstream of the pair). One refinement: attach iron to the U plane only,
  vacuum to the V plane, since U and V are 2.35 mm apart with no iron between
  them. Implement by checking `pi.plane`:
  ```cpp
  if (pi.plane == 0 && pi.layerInDet > 0) attach iron+scifi
  else if (pi.plane == 1)                 attach scifi only
  ```

### Step 3 — verify

Three checkpoints, each runnable in < 1 minute:

1. **Static check:** add an `info()` log printing each surface's
   `surface->surfaceMaterial()->materialSlab(Vector2{0,0}).thickness()` and
   `.material().X0()` for layer 0 of each subdetector. Eyeball values against
   PDG (Si X₀ ≈ 93.7 mm, Fe X₀ ≈ 17.6 mm, W X₀ ≈ 3.5 mm).

2. **Single-muon pipeline:** re-run
   `gaudi_jobs/muon_pipeline/muon_pipeline.sh` and the
   `diagnose_mtc_curvature.py` diagnostic. Expected outcomes:
   - sign-flip events drop from ~12% to ~0%.
   - median |R_truth / R_reco| drops from ~1.1 to within 1.0 ± 0.05.
   - reco track at MTC exit lies within tens of mm of truth (currently
     hundreds to thousands).
   - q/p prior is still loose so refine seed covariance in a follow-up.

3. **Compare against a no-MS reference:** temporarily set
   `multipleScattering=false, energyLoss=false` at
   `ACTSProtoTracker.cpp:1340-1341` and confirm the over-curvature returns —
   isolates whether the improvement is genuinely from MS or from a side
   effect of the geometry change.

### Risks of Path B

- Material is attached at the surface, not on a separate passive layer between
  SciFi planes. ACTS' MS noise contribution is the same to first order (it's
  a thickness × 1/X₀ integral, doesn't care about the geometric placement
  within the segment), but the location of the *energy loss* differs from
  what Geant4 simulated. Acceptable for tracking but should be noted in any
  future precision study.
- Tungsten / Silicon "averageLayers" combination loses the distinction
  between SiPad's absorber and sensor. For pure tracking it's adequate; for
  shower-shape analysis it isn't.
- Adding material can change CKF measurement-selector chi² distributions.
  May need to retune `Chi2CutOff` (`job4_tracking.py:101`).

### Rollback

Single git revert of `ACTSGeoSvc.cpp` restores prior behaviour. No other
files touched.

---

## Path A — full `Acts::convertDD4hepDetector` migration

### Why the bigger rework is worth doing eventually

Path B is a graft: it duplicates material knowledge (the X₀, ρ for Iron,
Tungsten, Silicon, Polystyrene) between the DD4hep compact XML and
`ACTSGeoSvc.cpp` even though we look the values up dynamically. As soon as
someone changes a layer thickness or adds a fourth subdetector, Path B drifts.

Path A makes the compact XML the single source of truth for both Geant4 *and*
ACTS by handing the geometry construction over to ACTS itself via the
DD4hep integration in `ActsExamples`.

### What `Acts::convertDD4hepDetector` does

(Reference implementation:
`Examples/Detectors/DD4hepDetector/src/DD4hepDetector.cpp` in the user's
private ACTS checkout, and the upstream
`Examples/Detectors/DD4hepDetector/include/.../DD4hepDetector.hpp`. The
canonical free function in 44.x is `ActsExamples::DD4hep::convertDD4hepDetector`;
its inputs are a `dd4hep::DetElement` (the world) and binning hints, its output
is a `std::shared_ptr<const Acts::TrackingGeometry>`.)

Under the hood it:

1. Walks the `DetElement` hierarchy starting from the world.
2. Consults each `DetElement`'s `Acts::ActsExtension` to decide whether the
   element is a `Layer`, a `Module`, a `Sensitive`, a `Passive`, a `Beampipe`,
   etc.
3. For each `Sensitive` `DetElement`, builds an `Acts::Surface` whose pose
   comes from `DetElement::nominal().worldTransformation()` and whose
   material comes from `DetElement::volume().material()` via
   `DD4hepDetectorElement`.
4. For each `Layer` `DetElement`, calls `Acts::DD4hepLayerBuilder` to collect
   that layer's sensitive surfaces and any passive material slabs declared
   inside it.
5. For each `Volume` `DetElement`, builds an `Acts::TrackingVolume` via
   `Acts::CuboidVolumeBuilder` or `CylinderVolumeBuilder` depending on shape.
6. Connects the volume tree into a `Acts::TrackingGeometry`.

Crucially, the **material on each surface is automatically pulled from the
DD4hep `Material` of the parent volume**, with no extra code. SiPad's
tungsten absorbers, if declared as their own `Passive` `DetElement`s,
become genuine passive `Layer`s in the tracking geometry with their own
`HomogeneousSurfaceMaterial`.

### Pre-requisite: annotate the plugins

`Acts::convertDD4hepDetector` needs each `DetElement` to carry an
`Acts::ActsExtension` declaring its role. None of our three plugins currently
do this. Each needs an audit and an update:

#### `detector_plugin/SiTargetDetector.cpp`

For every sensitive Si strip plane:
```cpp
auto* ext = new Acts::ActsExtension();
ext->addType("sensitive", "detector");
ext->addType("axes", "definitions", "XYZ");  // or whatever maps DD4hep→ACTS local
strip_plane_detElement.addExtension<Acts::ActsExtension>(ext);
```
For the SiTarget station envelope:
```cpp
auto* envExt = new Acts::ActsExtension();
envExt->addType("layer", "detector");
station_detElement.addExtension<Acts::ActsExtension>(envExt);
```

#### `detector_plugin/SiPadDetector.cpp`

Same pattern for sensors. Add a `passive`-tagged `DetElement` per tungsten
absorber slab:
```cpp
auto* tungstenExt = new Acts::ActsExtension();
tungstenExt->addType("passive", "detector");
absorber_detElement.addExtension<Acts::ActsExtension>(tungstenExt);
```
This is the step that makes the W absorbers visible to ACTS as MS-contributing
layers instead of being invisible Geant4-only volumes.

#### `detector_plugin/MTCDetector.cpp`

Sensitive: SciFi U/V planes and scintillator pads. Passive: each 50 mm iron
slab. Layer: each MTC layer (one SciFi U + one SciFi V + scintillator + iron).
Volume: each MTC station (MTC40, MTC50, MTC60).

The compact XML must also declare the iron material against each iron
volume's `<material name="Iron"/>` (already done for Geant4); no XML changes
needed if so.

### Replace `ACTSGeoSvc::initialize` body

After all three plugins are annotated:

```cpp
StatusCode ACTSGeoSvc::initialize() {
  // ... existing DD4hep load ...
  auto& desc = dd4hep::Detector::getInstance();

  // Build TrackingGeometry from the annotated DetElement tree.
  m_trackingGeometry = ActsExamples::DD4hep::convertDD4hepDetector(
      desc.world(),
      Acts::Logging::INFO,
      Acts::BinningType::arbitrary,   // r-binning (not used here, flat geom)
      Acts::BinningType::arbitrary,   // z-binning
      Acts::BinningType::arbitrary);

  // Cache surfaces by detector address for surfaceByAddress() — unchanged.
  // ...
  return StatusCode::SUCCESS;
}
```

Drop the `extractPlanes` lambda (~250 LoC), the `PlaneInfo` struct, the
`SNDDetectorElement` class, and the explicit `PlaneLayer::create` /
`LayerArrayCreator` / `TrackingVolume` block. About 350 LoC of custom geometry
construction becomes a single function call.

### Caveats

1. **The `CuboidVolumeBuilder` offset bug.** The current code bypasses it for
   the reason noted at `ACTSGeoSvc.cpp:240-243`. `convertDD4hepDetector` uses
   the same `CuboidVolumeBuilder`. Either:
   - re-verify against Acts 44.3 that the bug is fixed, or
   - submit a patch upstream first, or
   - keep Path B if the bug recurs.

2. **The `SNDFixedNavigator` workaround.** This was needed because the
   standard `Acts::Navigator` failed to traverse our single-volume layout
   (see `docs/known_issues.md` "Root cause 3"). With a proper multi-volume
   hierarchy out of `convertDD4hepDetector`, standard `Acts::Navigator`
   should work and `SNDFixedNavigator` can be retired. Confirm with a test
   propagation before deleting it.

3. **`surfaceByAddress(detID, station, layer, plane)`** in `ISNDGeoSvc` is
   currently keyed by our own `PlaneInfo` indexing. After migration the
   address-keying must come from DD4hep cellIDs / `DetElement` paths. Add a
   compatibility shim that maps the old `(detID, station, layer, plane)`
   tuple to a `dd4hep::DetElement` (or a `GeometryIdentifier`) — otherwise
   `MTCSciFiMeasConverter.cpp:184` and `ACTSProtoTracker.cpp:1008` break.

4. **Build system.** `Acts::convertDD4hepDetector` lives in
   `ActsExamples::DD4hep`, which is an Examples-tier library — usually not
   installed by Acts' default release. We may need to either:
   - link to the user's private Acts checkout
     (`/afs/cern.ch/user/e/edursov/private/SIMULATIONS/my_FairShip/...`), or
   - pull the relevant code into our tree as `gaudi_source/SNDDD4hepConverter.cpp`
     (~500 LoC, MPL-2 licensed), or
   - upstream a request to expose this from the Acts core release in
     `2026-02-01`+.

### Verification (Path A)

Run the diagnostic and compare against Path B baseline. Required pass
criteria:
- All Path B verification criteria still met.
- Material X₀ per surface read directly from DD4hep matches what Path B
  hard-coded (use the static dump in Path B Step 3 ① as the reference).
- `Acts::TrackingGeometry::visitSurfaces` enumerates all ~150 surfaces with
  non-null `surfaceMaterial()`.
- ddsim + Gaudi pipeline produces identical hit counts pre- and post-migration
  (the geometry seen by Geant4 is unchanged; only the ACTS view changed).

### Rollback

A flag-based switch in `ACTSGeoSvc` lets Path A and Path B coexist during
transition:
```cpp
Gaudi::Property<bool> m_useConverter{this, "UseDD4hepConverter", false,
    "If true, use Acts::convertDD4hepDetector instead of manual geometry build."};
```
Default `false` (= Path B / current). Flip to `true` only after the migration
is verified end-to-end. Removes the rollback risk during the longer migration
work.

---

## Order of work

1. Implement Path B Step 1-3, verify with diagnostic. **Closes the
   over-curvature bug.**
2. (Optional, after step 1 lands) tighten `seedCov(eBoundQOverP)` to ~0.04 and
   tune `Chi2CutOff`.
3. Audit plugins for `DetElement` granularity. Document gaps.
4. Add `ActsExtension` annotations to the three plugins.
5. Add `m_useConverter` switch to `ACTSGeoSvc` plus the
   `convertDD4hepDetector` code path.
6. Verify against Path B baseline.
7. Retire `SNDFixedNavigator`, `SNDDetectorElement`, and the manual TGeo walk.
8. Update `docs/known_issues.md` to mark the "no surface material" risk as
   fixed and remove the latent-risk note for the sign convention if the
   Convert path obviates it.

---

## Open questions

- **DD4hep `Material` API units.** ✅ Confirmed: `radLength()` and `intLength()`
  return cm in the 2026-02-01 stack. The `* 10.0` factor in `makeSlab` is correct.
- **MTC iron between U/V planes.** ✅ Confirmed: no iron between U and V.
  `MTCDetector.xml` layer order: Iron(50mm) → Iron(3mm) → Si-U → Air(1mm) → Si-V → Iron(3mm) → Scint.
  Iron(53mm total) attaches to U plane only; V plane gets scintillator only.
- **SiPad tungsten thickness.** ✅ Confirmed 3.5mm from `SND_compact.xml:43`
  (`Ecal_WThickness = 3.5*mm`). Material name: `TungstenDens1910`.
- **`MaterialSlab::averageLayers` availability.** ✅ Does not exist in 44.3.
  Use `Acts::MaterialSlab::combineLayers(slabA, slabB)` — confirmed in header.
- **MTC SciFi material name.** ✅ The compact XML uses `material="Silicon"` for the
  SciFi planes (geometry placeholder), but `"Scintillator"` is used in the
  material attachment to reflect the physically correct fibre material.
- **CuboidVolumeBuilder offset bug.** Still open — relevant only for Path A.

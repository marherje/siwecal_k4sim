# MTC Over-Curvature Investigation

**Status: RESOLVED 2026-05-17** — Path B material implementation in `gaudi_source/ACTSGeoSvc.cpp`.
See resolution summary below and `docs/known_issues.md` "No surface material" entry.

**Pipeline:** `gaudi_jobs/muon_pipeline/` with input
`output/output_mu-_xyz_5.0_-82.5_-1000.0_dir_0_0_1_E5_100_events.edm4hep.root`
(50 mu- events, 5 GeV, gun at (5, -82.5, -1000) mm, direction (0,0,1)).

**Diagnostic script:** [`diagnose_mtc_curvature.py`](diagnose_mtc_curvature.py)
Plots in [`plots_curvature/`](plots_curvature/).

---

## Resolution (2026-05-17)

**Root cause confirmed:** `ACTSGeoSvc` attached no `ISurfaceMaterial` to any surface.
The CKF MS noise term reads `surface->surfaceMaterial()` and contributes zero when null.
With MS disabled (effectively), each 53 mm of MTC iron was treated as vacuum and all MS
kicks were absorbed into the fitted q/p, over-tightening the reconstructed radius by 3–10×.

**Fix:** Path B from `docs/acts_material_migration.md` — `makeSlab` helper added to
`ACTSGeoSvc.cpp` (anonymous namespace) queries `dd4hep::Detector::material(name)` for
X₀, L₀, ρ and attaches `HomogeneousSurfaceMaterial` to each surface after construction.
Iron(53mm) + Scintillator(1.35mm) on MTC U planes; plain Scintillator on V planes;
TungstenDens1910(3.5mm) + Silicon on SiTarget/SiPad planes.

**Post-fix results (same 50-event input):**

| Metric | Before | After |
|--------|--------|-------|
| Median \|R_truth / R_reco\| | 3–10× | **0.94** |
| Median reco \|dx\| through MTC | 1600–3400 mm | **1–6 mm** |
| Median \|dx\|_max through MTC | ~1700–3400 mm | **10–20 mm** |
| Sign agreement | 21/24 = 87.5% | 22/26 = 85% |

The 4 sign-flip outliers in the post-fix run are seeding ghosts (tracks sharing
hits with the primary muon); they are not material-related and were present before.

---

## Symptom

For a 5 GeV mu- in 1.7 T MTC iron, the truth and reco trajectories agree at the
MTC entrance but diverge dramatically toward the MTC exit. The reco track keeps
bending in the same direction as truth (mostly), but with much stronger curvature.

| Window | Truth Δx over MTC | Reco Δx over MTC | Reco / Truth |
| ------:| -----------------:| ----------------:| ------------:|
|   0    |   +262 mm         |   +1642 mm       |  6.3×        |
|   1    |   +265 mm         |   +1686 mm       |  6.4×        |
|   2    |   +201 mm         |   −1879 mm       |  −9.3×       |
|   4    |   +257 mm         |   +2517 mm       |  9.8×        |
|   5    |   +255 mm         |   +1621 mm       |  6.4×        |
|   6    |   +202 mm         |   +1998 mm       |  9.9×        |
|   7    |   +243 mm         |   +3394 mm       |  14×         |
|   8    |   +205 mm         |   +1570 mm       |  7.7×        |
|   9    |   +256 mm         |   +2615 mm       |  10×         |

Truth Δx ≈ 250 mm matches the analytic expectation for a 5 GeV mu-
(R = p/(0.3·B) ≈ 9.8 m; ≈ 76 mrad bending per station; ≈ 250 mm net deflection
over the three stations).
Effective reco radius from the deflection: R_reco ≈ 1.4–2.9 m, giving
p_reco ≈ 0.5–1.5 GeV — **3–10× too low**. Sign is correct in 21/24 fitted tracks
(87.5 %); the remaining 12.5 % bend in the opposite direction (window 2 above).

## Where the divergence happens

`diagnose_mtc_curvature.py` shows the reco track entering MTC at the truth
position (within 1–4 mm) and then progressively peeling away from the truth
through the three stations. The over-curvature is **cumulative**, not a single
bad seed: the KF starts correctly, then drifts.

```
                                    MTC entrance   MTC exit
truth                     5 mm  →  −3..14 mm   →   210..270 mm
reco (sorted by chi²)     5 mm  →    1..12 mm   →  1600..3400 mm  (or −1900 mm)
```

## Suspect components (ranked by likelihood)

### 1. Wrong charge sign on the seed q/p — `ACTSProtoTracker.cpp:1207-1208`

```cpp
const double seedQoverP =
    1.0 / (m_seedMomentum.value() * Acts::UnitConstants::GeV);   // always +1/p
…
Acts::BoundTrackParameters::create(
    gctx, sfSeed->getSharedPtr(), seedPos4, seedDir, seedQoverP,
    seedCov, Acts::ParticleHypothesis::muon());
```

The truth particle is a **mu-** (charge −1, q/p < 0), but the seed q/p is
**positive**. `Acts::ParticleHypothesis::muon()` carries only the absolute
charge (1·e); the sign comes from the qOverP scalar. With the (loose) qOverP
prior variance of 10 (1/GeV²) the CKF can numerically flip the sign in most
cases (21/24 events here) but it cannot do so smoothly: q/p → 0 corresponds to
p → ∞, so the optimiser can settle on a far-from-truth |q/p| even after a sign
flip. In the 3/24 events where the flip fails the track bends the wrong way
entirely. In the rest the flipped value lands at |q/p| larger than truth, i.e.
**too small a momentum → too tight a radius → over-curvature**.

### 2. MultipleScattering and EnergyLoss disabled — `ACTSProtoTracker.cpp:1336-1341`

```cpp
Acts::CombinatorialKalmanFilterOptions<SNDTrackContainer> ckfOptions(
    gctx, m_mctx, std::cref(m_cctx),
    ckfExtensions, pOptions,
    /*multipleScattering=*/false,
    /*energyLoss=*/false);
```

For a 5 GeV muon traversing 750 mm of iron per station (≈ 42.6 X₀) the multiple-
scattering angle θ_MS ≈ 13.6 MeV / p · √(x/X₀) · (1 + 0.038·ln(x/X₀)) ≈ 18 mrad
per station, vs ≈ 76 mrad of true bending. The KF without an MS process-noise
term treats each MS kick as a real direction change and "absorbs" it into the
fitted q/p. Over 45 iron slabs (3 stations × 15 layers) this drives the fitted
momentum systematically downward (curvature upward).

### 3. Seed momentum 10 GeV with truth 5 GeV — `job4_tracking.py:97`

```python
proto.SeedMomentum = 10.0   # GeV
```

A 2× wrong seed is normally not a problem (qOverP variance is huge → KF
re-estimates from measurements). It is, however, the wrong *sign* of correction
to debug: with seedMomentum > truth and KF correctly updating, p_fit should
approach 5 GeV from above. Instead p_fit is **below** truth, which again points
to either a sign issue (1) or a noise-model issue (2) — not the seed-p value
itself.

### 4. (Not the cause) B-field strength / direction

Verified consistent between Geant4 (`detector_plugin/MTCDetector.cpp`,
`MTCIronField::fieldComponents`) and ACTS (`ACTSProtoTracker.cpp`,
`IronSlabBField::getField`). Both apply +1.7 T along the global Y axis inside
the 50 mm outer-iron slabs only. ACTS converts via `Acts::UnitConstants::T`
correctly. Geometry slab extents agree to ~1 mm.

## Suggested fixes (in order)

1. **Sign-correct the seed q/p**: pass a negative q/p for negatively charged
   hypotheses. Easiest patch:
   ```cpp
   const int    seedSign  = m_seedCharge.value();  // new Gaudi property, default −1 (mu-)
   const double seedQoverP = seedSign /
       (m_seedMomentum.value() * Acts::UnitConstants::GeV);
   ```
   Re-run the muon pipeline and re-evaluate. Expect the sign-flip-failure events
   (windows 2, …) to disappear immediately, and the |Δx|_exit to drop.
2. **Enable MS modelling** (and optionally energy loss) in the CKF:
   ```cpp
   /*multipleScattering=*/true,
   /*energyLoss=*/true);
   ```
   This gives the KF the correct process noise so it stops over-absorbing
   scattering kicks into q/p.
3. **Lower SeedMomentum to a realistic prior** (5 GeV for this gun). Less
   critical than (1)/(2) but reduces the dynamic range the KF has to traverse.

## How to reproduce

```bash
source init_key4ship.sh
cd gaudi_jobs/muon_pipeline
# regenerate ShipHits.root if needed
source muon_pipeline.sh
# run the diagnostic
python diagnose_mtc_curvature.py --max-windows 10
# plots: plots_curvature/mtc_curvature_w*.pdf
```

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
          (DIGI_MODE=real|both swaps in / adds RealDigitizer -- see Job 3)
  → job4: ShowerTagger + SiPadMeasConverter + ACTSProtoTracker
                              → digitized.edm4hep.root  (ACTSTracks + EMShowers
                                 + SiPadMeasurements, written back into the
                                 SAME file via a temp-output + swap in the
                                 pipeline script -- see "Single-file output"
                                 below)
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
`RealDigitizer.cpp`, `include/CellShaping.hh`, `DetectorFlipper.cpp`,
`ChannelMapper.cpp`
**Config:** `gaudi_jobs/pid2026_common/job3_digitize.py` (shared)

There are two digitisation chains. `DIGI_MODE` picks which one(s) run. **Since
2026-09-21 the default is `real`**, the chain anchored to the test beam; the
historical `simple` chain has to be asked for. In the same change
`HIT_SELECTION` defaults to `chip` and `CALIB_THRESHOLD` lost its default:
job3 refuses to run without it, because every calibrated number depends on the
threshold set and nothing may guess it (the launchers export it, a one-off run
states it on the command line).

```
DIGI_MODE=real    (default)
SiPadHits → RealDigitizer     → SiPadHitsRealDigi  <-- tracking input (pre-flip)
                              + SiPadHitsRealDigitizedEnergy [MIP]
                              + SiPadHitsRealDigitizedTime   [ns]
                              + SiPadHitsRealDigitizedTrigger (HIT_SELECTION=chip)
          → DetectorFlipper   → SiPadHitsRealFlipped
          → ChannelMapper     → SiPadHitsRealMapped + SiPadHitsRealMasked
          → AdcDigitizer      → SiPadHitsRealAdc + ...AdcHigh/Low/Kept   <-- ecal tree

DIGI_MODE=simple
SiPadHits → GeV2MIPConversion → SiPadHitsMIP
          → BasicDigitizer    → SiPadHitsDigi      <-- tracking input (pre-flip)
          → DetectorFlipper   → SiPadHitsFlipped
          → ChannelMapper     → SiPadHitsMapped + SiPadHitsMasked

DIGI_MODE=both    both of the above, same input, same output file
```

The launchers follow the mode: `generic_condor_beam_chunk.sh`,
`reprocess_chunk.sh` and `1_mu_pipeline.sh` track on `SiPadHitsRealDigi` and
convert `SiPadHitsRealAdc` (with `SiPadHitsRealMasked`) under `real`, and the
`SiPadHitsDigi` / `SiPadHitsMapped` pair under `simple`; `run_pid_sim.sh`
defaults to `--collection SiPadHitsRealAdc`. The geometry defaults moved with
it: `generic_condor_beam.sh` and the chunk launcher simulate
`SND_compact_beamline.xml` with the beam at `BEAM_Z_MM=-2650` and tag the
output `_beamline_real`; the bare detector is
`COMPACT_FILE=SND_compact.xml BEAM_Z_MM=-2000 OUTPUT_TAG=`.

| Algorithm | Role |
|---|---|
| `GeV2MIPConversion` | Energy in GeV → MIPs. `MIPValues` takes one value per layer (from `mip_extraction_pipeline`); `MIPValue` is the scalar fallback |
| `BasicDigitizer` | Applies the MIP `Threshold`, dropping hits below it |
| `RealDigitizer` | `DigitizationMode='real'`: per-cell CR-RC shaping over the `CaloHitContributions`. The fast channel gives the trigger time at its first threshold crossing, the slow channel is sampled at `triggerTime + DelayNs`. Takes the same per-layer `MIPValues` as `GeV2MIPConversion`, and its trigger model (`Threshold`, `FastNoiseMIP`, `TriggerEfficiency`) from the threshold set — see "One threshold input" below. `DigitizationMode='simple'` is the `BasicDigitizer` threshold cut |
| `DetectorFlipper` | Rewrites hit z into the **test-beam frame** using `mappings/slab_z_positions.yml` — the single source of truth for the per-slab z, shared with the event viewer |
| `ChannelMapper` | Cell IDs → test-beam format (`system:8,slab:8,chip:16,channel:8,sca:8`) via the pad maps in `mappings/`, and masks dead channels from the muon calibration tree (`CalibThreshold`, e.g. `th230`, masks ~3.5% of channels) |

**Tracking reads `SiPadHitsDigi`, before `DetectorFlipper`**: the flip moves the
hit z into the test-beam frame, which no longer matches the ACTS surfaces built
from the compact XML.

### One threshold input — `CALIB_THRESHOLD`

Everything that depends on the trigger threshold is driven by **one** environment
variable, `CALIB_THRESHOLD` (`th210` / `th220` / `th230`, no default — required):

```bash
DIGI_MODE=both CALIB_THRESHOLD=th220 INPUT_FILE=... k4run gaudi_jobs/pid2026_common/job3_digitize.py
```

- `ChannelMapper.CalibThreshold` — which `mips/<th>/` table masks the dead channels
- `RealDigitizer.Threshold` / `.FastNoiseMIP` / `.TriggerEfficiency` — the trigger
  model of that set, read from `mappings/trigger_turnon.yml`

The trigger model is **measured**, not assumed. `hitbit_high` in the decoded
chunks is the fast shaper's discriminator, and `EventBuilder::bestScaPerChannel`
drops any channel whose bit never fired, so its turn-on IS the data's hit
selection. `analysis/trigger_turnon.py` fits it per threshold set with
`P(x) = plateau · ½(1 + erf((x − mu)/(√2 σ)))`, whose three parameters map onto
the three properties one for one:

| set | run | threshold [ADC] | σ [ADC] | threshold [MIP] | σ [MIP] | plateau |
|---|---|---|---|---|---|---|
| th210 | eudaq 253 | 16.26 | 2.64 | 0.886 | 0.150 | 0.964 |
| th220 | run 72 | 19.67 | 4.05 | 0.872 | 0.165 | 0.852 |
| th230 | run 13 | 25.36 | 4.66 | 0.781 | 0.141 | 0.862 |

Before this the digitiser applied a hardcoded 0.5 MIP step with 100% efficiency,
whatever the sample was meant to represent.

**The ADC column is the trustworthy one.** The discriminator level rises
monotonically with the DAC (16.3 → 19.7 → 25.4 ADC) while the pedestals and the
LG→HG anchor stay put, which is the calibration-free statement that the DAC moves
the threshold and nothing else. The MIP column divides each by *that set's own MPV
table*, so it inherits whatever bias the table has — which is exactly why th230
comes out *lowest* in MIP despite having the highest discriminator. When the ADC
model lands (pedestal, gain, quantisation), the digitiser should take
`threshold_adc` and stop going through MIP at all.

Two implementation notes, both of which were bugs the turn-on exposed:

- The pre-cut that skips the shaping for small cells now keeps a 5σ margin below
  the threshold. Cutting exactly at the threshold removed every cell below it
  *before* the noise was added, so the simulated turn-on could only ever be a
  step, whatever `FastNoiseMIP` said.
- The threshold decision and the trigger-time search now run on the **same**
  waveform. The noise used to be added only to the peak while the crossing was
  hunted on the clean curve, so a cell whose noisy peak passed but whose clean
  peak did not never found a crossing and was dropped anyway.

### The ADC model — `AdcDigitizer`

**File:** `gaudi_source/AdcDigitizer.cpp`, after `ChannelMapper` (which is what
knows the test-beam channel, and therefore which pedestal and MPV apply).
`ADC_MODEL=0` turns it off and leaves the chain in raw digitised MIP.

Until this existed the simulated chain worked in MIP from end to end, so it never
had the thing that stands between charge and energy in the data: **a
calibration**. The test beam's energy is `(adc_high − pedestal) / MPV`, with the
MPV measured from a muon run. Dividing the simulation by its own clean MIP scale
and the data by a measured one does not compare the same quantity.

So the algorithm does what the detector *and* the reconstruction do, in order:

1. **MIP → ADC** with `AdcPerMip`, the absolute charge scale, **fitted per
   threshold set against that set's own test-beam data** (see "Calibrating the
   ADC scale" below). Only the channel-to-channel *variation* comes from a table
   (`GainShapeThreshold`, normalised to its median): the scale is data.
2. **The readout**: add the pedestal, round to integer ADC, clamp the high gain at
   `AdcHighMax` — **measured** from the upper edge of that set's `hit_hg`
   spectrum (2046 / 2144 / 2062 ADC for th210 / th220 / th230; the preamp rolls
   over rather than clipping, so a hard ceiling is an approximation) — and carry
   the low gain through the anchor line `adc_low = k·adc_high + c`.
3. **ADC → MIP** with the **reconstruction** table (`CalibThreshold`), low-gain
   recovery above `SaturationAdc` included — a port of `EventBuilder::buildHit`.

The output hits (`SiPadHitsRealAdc`) carry the *reconstructed* MIP, and the
parallel `SiPadHitsRealAdcHigh` / `...AdcLow` carry the pedestal-subtracted ADC,
which `sim_to_ecal_tree.py` puts in `hit_hg` / `hit_lg`. The simulated ecal tree
is the same object as the data's, ADC included. The algorithm is strictly 1:1 so
the parallel masking / time / fast-peak collections stay index-aligned.

Every per-threshold number lives in **`mappings/digi_calibration.yml`** — the
trigger (`threshold_adc`, `sigma_mip`, `efficiency`), the scale (`adc_per_mip`),
the ceiling (`adc_high_max`) — each with its provenance, and `CALIB_THRESHOLD`
picks the entry. The trigger threshold `RealDigitizer` applies is
`threshold_adc / adc_per_mip`, so it moves with the scale. The trigger and the
dynamic range are that set's own; the gain is **one number for the three sets**
(19.5 ADC/MIP since 22 Sep, see below) — the preamplifier does not know the
trigger DAC. No MIP-table ratio is inherited between sets
(`test_every_threshold_set_is_self_contained` fails if one comes back).

### Calibrating the ADC scale — `analysis/fit_adc_scale.py`

Per threshold set, three numbers, each from the data that measures it and none
inherited from another set:

1. **The gain** (`adc_per_mip`) from a **muon run**: the `adc_per_mip` whose
   simulated per-hit high-gain spectrum has the data's shape
   (`analysis/muon_gain_scan.py`), cross-checked by the tag-and-probe layer
   efficiency of muon tracks. The MIP is the gain by definition, with no shower
   physics in between — **but only where the discriminator sits below the MIP
   peak** (th210). Above it (th220, th230) the surviving Landau tail is
   scale-free and the scan measures nothing (see *One gain*, below); the one
   value measured at th210 serves every set.
2. **The dynamic range** (`adc_high_max`, `saturation_order`) from an **electron
   run** of that set, `analysis/compare_gains.py --write`: the low gain is linear
   over the whole range, so `(hit_lg − c)/k` is the true charge and `hit_hg`
   against it is the high gain's response curve, fitted with
   `S(q) = q / (1 + (q/A)^n)^(1/n)`. Only electrons reach the knee. This fit
   involves no simulation at all. The same pass measures how the data's low
   gain scatters about the anchor line (`lg_noise_adc`, 3.5–4.0 LG-ADC of
   pedestal width, and `lg_gain_spread`, ~7 % growing with the amplitude — the
   channel-to-channel dispersion of the gain ratio) and `AdcDigitizer` draws
   them per hit (`LowGainNoiseAdc`, `LowGainSpread`): the simulated
   `hit_lg`-vs-`hit_hg` band then has the data's width. Both are symmetric, so
   they widen the energy of the hits reconstructed from the low gain (above
   `SaturationAdc`) and move nothing on average — the 52 GeV event sum is the
   same to four digits with and without them.
3. **Every electron energy is then a validation**, nothing refitted.

A shower anchor for the gain — tried first — is *not* a gain measurement: it
absorbs whatever the simulation gets wrong about the shower. Anchoring th220 on
74 GeV showers gave 15.6 ADC/MIP and failed the muon check outright (5 simulated
hits per muon event against 20: the threshold landed above the MIP peak).

The beam positions differ per set — th230 ≈ (−42, +51), th220 ≈ (−40, +37),
th210 ≈ (+24, −51) mm for electrons, and the muon runs sit elsewhere again —
measured from each run's `hit_hg`-weighted barycentre, so the productions are
position-matched (`launch_beam.sh` / `generic_condor_beam.sh`). The muon runs of
th230 (run 4) and th210 (eudaq 151) had no ecal tree and were reconstructed with
`siwecal-tb2026/gaudi_jobs/condor/generate_reco_dag.py`.

Two things the muon anchor turned up, both now in the code:

- **`RealDigitizer` applied a second threshold on the slow sample** after the
  fast discriminator had fired. The chip does not: `hitbit_high` is the only
  hit selection, and `bestScaPerChannel` asks for nothing else. The product of
  the two cuts put the simulated 50% point at 1.0 MIP instead of
  `threshold_adc / adc_per_mip` and narrowed the S-curve to two thirds of its
  width. Removed — the simulated turn-on now lands exactly where it is
  configured (0.85 / 0.91 / 1.16 MIP for th210 / th220 / th230, width 0.15–0.17).
- **`EcalEventBuilder` keeps an event only if `MinSlabsHit = 10` slabs fired.**
  Irrelevant for a shower, decisive for a muon: it is a selection on the Landau
  fluctuations of its fifteen crossings. The anchor applies it to the simulated
  sample too (`--min-layers 10`), plus `--max-hits 30` against pile-up (43% of
  run 4's events are above 30 hits and they move the peak by 2 ADC).

At th230 the discriminator (25.4 ADC) sits *above* the MIP, so what survives is
the Landau tail — which is scale-free — and the mode of the surviving spectrum
is not monotonic in the gain (26.6 / 27.6 / 27.5 / 28.4 / 26.6 ADC for gain
20 / 21 / 22 / 23 / 24): a fixed-point fit on it lands wherever it starts. The
gain is therefore taken from the **shape** of the whole spectrum just above the
threshold (`analysis/muon_gain_scan.py`: χ² of the unit-normalised spectra over
[threshold + 1, 80) ADC against a 1 ADC/MIP scan), which has one minimum. The
same scan is run on the other two sets as a cross-check of the mode fit.

Measured (2026-09-18), 10 000-event muon samples, 1000-event electron samples;
the electron rows are validations with nothing refitted:

| set | gain [ADC/MIP], scan (mode fit) | muon run | A [ADC] | n | discriminator [MIP] | electron point | sim/data ADC peak | sim/data hits |
|---|---|---|---|---|---|---|---|---|
| th210 | **19.0** ± 1.4 (18.3) | eudaq 151 | 1951 ± 26 | 4.9 | 0.85 | e⁻ 52 GeV, eudaq 287 | **1.07** | 1.13 |
| th210 | | | | | | e⁻ 74 GeV, eudaq 286 | **1.15** | 1.45 |
| th220 | **21.5** ± 1.4 (21.1) | run 85 | 1918 ± 16 | 6.1 | 0.91 | e⁻ 74 GeV, run 72 | **1.33** | 1.40 |
| th230 | **21.9** ± 1.5 (22.0) | run 4 | 1930 ± 41 | 5.8 | 1.16 | e⁻ 20 GeV, run 20 | **1.12** | 1.05 |
| th230 | | | | | | e⁻ 52 GeV, run 13 | **1.12** | 1.12 |

(the ± is where the scan's χ² doubles, not a 1σ; the ADC peak is the mode of
the per-event `sum_hg`, the hits are `nhit_chan` means.)

Reading:

- **The gain does not follow the DAC**: 19.0 / 21.5 / 21.9 ADC/MIP, the same
  within the scan's resolution, as it must be for a preamplifier setting the
  DAC does not touch. th210 is the lowest and is also the eudaq run at a
  different position, so a few percent of channel-to-channel gain is plausible.
- **The th230 MIP table (median 32.4 ADC/MIP) is 1.48× the gain the muons
  give.** With 21.9 ADC/MIP the th230 discriminator is at **1.16 MIP**, above
  the peak — consistent with the table having been fitted on a truncated
  Landau, which is where the ×1.47 first showed up.
- Validations: **within 7–15% at 20 and 52 GeV** at th230 and th210, and 1.33
  at 74 GeV for th220 (1.15 for th210). The excess grows with energy at fixed
  gain and roll-over, so it is not the ADC scale nor the high-gain ceiling: the
  data's response per GeV falls with energy and the simulation's does not (the
  occupancy loss already seen in MIP/GeV: 34 → 30 → 23 from 52 to 99 GeV).
- The simulated hit count is high by 5–13% at 20–52 GeV and 40–45% at 74 GeV,
  the same trend.
#### One gain for the three sets (22 Sep) — supersedes the per-set gains above

The per-set gains 19.0 / 21.5 / 21.9 (19.07 / 21.33 / 22.39 after the chip
selection) were not measurements above the MIP. `muon_gain_scan.py` used a
Pearson χ² with the data as the only expectation; that statistic has an
expected value of n_bins / N_sim_hits even for identical shapes, and N_sim_hits
grows with the gain because more simulated events pass `MinSlabsHit ≥ 10`. At
th230 the "minimum at 22.4" was that floor — the two curves coincide for gains
19–23 and only ≥ 24 is excluded (`final_v3/plots/gain_scan_floor.png`,
`gain_scan_check.py`). At th220 the minimum at 21.3 tracks the 0.85 plateau
measured on showers and applied to muons; with th210's 0.96 it would sit at
~19. The script now carries both Poisson variances and prints the floor.

Three estimates where the peak is visible: th210 spectrum shape with the
corrected χ² 18.9–19.1; tag-and-probe layer efficiency of muon tracks (data
0.867 at th210) 19.7; the th230 table's ×1.47 bias needs a threshold at
1.35 MPV (`mip_threshold_bias.py`) → 25.36 / 1.35 = 18.8. Adopted
**19.5 ± 0.5 ADC/MIP** for every set (`adc_per_mip` in
`digi_calibration.yml`, the per-set values kept as
`adc_per_mip_per_set_2026_09_18`). The discriminators are then 0.83 / 1.01 /
1.30 MIP (th210 / th220 / th230), slab 12 at 1.24 / 1.54 / 1.58, and the th230
table sits ×1.66 above the gain.

Result, same trees otherwise (basic beampipe, per-slab digitizer, chip / adc
selection), Gaussian core μ of the event ADC sum, hits > 0.5 MIP
(`Processed/adc_vs_tb/final_v3/`):

| set | point | sim/data ADC, per-set gain | **sim/data ADC, one gain** | hits | depth data / sim |
|---|---|---|---|---|---|
| th230 | e⁻ 20 GeV, run 20 | 1.16 | **1.00** | 0.98 | 5.49 / 5.75 |
| th230 | e⁻ 52 GeV, run 13 | 1.16 | **1.01** | 1.03 | 6.16 / 6.37 |
| th220 | e⁻ 74 GeV, run 72 | 1.26 | **1.16** | 1.13 | 6.71 / 6.59 |
| th210 | e⁻ 52 GeV, eudaq 287 | 1.08 | **1.10** | 1.05 | 6.23 / 6.38 |
| th210 | e⁻ 74 GeV, eudaq 286 | 1.13 | **1.16** | 1.10 | 6.80 / 6.66 |

The two 74 GeV sets now agree (the v2 spread was the two gains), and the
"excess of hits > 300 ADC in the back layers" is gone: it was 15 % on the ADC
axis acting on a steeply falling spectrum. What is left is a profile shift —
sim/data 0.72 at layer 0 rising to 1.0–1.2 at layers 9–14 with the total at
1.0, the simulated shower 0.25–0.3 layer deeper — which points at ~0.1–0.15 X₀
more in front of the box (XCET gas?) and/or ~5 % more X₀ per sampling cell, and
a rise with energy (1.00 → 1.16 at 74 GeV) that is the data-side loss at high
rate. Physics list (FTFP_BERT_EMZ, QGSP_BERT_EMZ) and a 0.05 mm range cut
change the total by 1–2 % and the profile not at all. Largest remaining
inconsistency between sets: 1.01 (th230) against 1.10 (th210) at 52 GeV — the
plateau, 0.96 in eudaq 253 and 0.82 in eudaq 287 at the same DAC, is measured
on showers and applied to muons; a tag-and-probe measurement on the muon runs
per set is the next step.

Earlier per-set reading, kept for the record:

- Muons: after the ≥ 10-slab selection the data carry 1.27 / 1.49 / 1.76 hits
  per fired slab (th210 / th220 / th230) against 1.17 / 1.20 / 1.44 in the
  simulation — extra neighbour or noise hits in the data that grow with the
  threshold, which is not what noise does. Open.

The three dynamic ranges agree within errors, as they should for a property of
the preamplifier rather than of the threshold. The DAC law from the three
turn-ons: **0.52 ADC per DAC unit** (run 20 for th230), residuals under 1 ADC.

All figures, the digitised trees, the gain-scan samples and this table are at
`/eos/experiment/drdcalo/siw-ecal/TB2026-06/Simulation/Processed/adc_vs_tb/final/`.

### Per-slab technology: slab 12 (the FEV11 chip-on-board) and layer 14 (650 µm)

With the same hit selection on both sides the per-layer sim/data ratio was flat
except at layer 12 (2.2) and layer 14 (1.4–1.6). Both are per-slab hardware facts
that nothing described:

- **Slab 12 is the FEV11 chip-on-board.** It ran at its own threshold DAC (243
  against the sets' 215–230; run book), which a per-slab S-curve confirms
  directly: discriminator at **30.8 / 30.1 / 24.2 ADC** for th230 / th220 / th210
  where the sets sit at 25.3 / 19.7 / 16.2 (`trigger_turnon.py --per-slab`). The
  truncated MIP spectrum inflates its table MPV (1.91 / 1.63 / 1.22 × the median
  in th210 / th220 / th230 — a real gain would not depend on the DAC), and the
  digitiser was inheriting that as *gain shape*.
- **Layer 14 has the 650 µm sensor** (geometry), but job3 normalised it with the
  500 µm MIP value (×1.34 in MIP) and the shape table added another ×1.25.

What changed:

| where | what |
|---|---|
| `mappings/slab_z_positions.yml` (both repos, identical) | technology block: `technologies`, `slab_technology` (slab 12 = `FEV11_COB`), `sensor_thickness_um` (650 at 14), `threshold_dac` (243 at 12). Loader `analysis/slab_description.py`. |
| `gaudi_jobs/mip_extraction_pipeline/` | **step 0**: `mip_extraction_pipeline.sh <muon ddsim>` writes `mappings/mip_values_sim.yml` (per-layer MIP peak, Landau⊗Gauss; layer 14 = 1.34 × the others); job3 refuses to run without it. The simple chain uses the per-layer values; the real chain uses the **charge unit** (the 500 µm value for every layer) because the ADC follows the charge and the reconstruction table turns it back into the data's MIP. |
| `ChannelMapper` | `PadMapSlabOverrides` (`"slab:path"`, from the description); `PadMapFileSlab12` kept as a deprecated alias. |
| `RealDigitizer` | `ThresholdPerLayer` / `FastNoiseMIPPerLayer` / `TriggerEfficiencyPerLayer`; job3 fills them from `slab_overrides` in `digi_calibration.yml` (in ADC, measured by `trigger_turnon.py --per-slab --write-config`). |
| `AdcDigitizer` | `GainShapeNormalisation = per-slab` (job3 default): each slab's table block divided by its own median, so a slab-level table offset never enters as gain. |
| geometry | layer 12 is its own `<layer>` block with `Ecal_WaferThickness_L12` / `Ecal_w_slab_gap_L12` (500 µm; one constant to switch); `analysis/tests/test_sensor_thickness.py` pins the XML against the YAML. |
| `siwecal-tb2026` | `SlabGeometry::fromYamlFile` reads the block (and no longer appends unknown lists to `w_thickness_mm`); `load_slab_technology()` in `siwecal_eventbuilder/geometry.py`; the `12:` pad-map override is derived from the YAML in `run_event_builder.py`, `run_full_pipeline_batch.py`, `generate_reco_dag.py`, the viewer and `add_xy_branches.py`; `PedestalMipCalibrator` has per-slab MIP windows (`SlabZFile`, `MipWindowSlabOverrides`). |

Result on the run-4 muons (per-hit peak, sim/data per layer): layer 0 1.27 → 1.04,
layer 12 1.24 → 0.92, layer 14 1.23 → 1.00; every layer within ±10 % except 6
(the data's own outlier) and 13. On the 52 GeV th230 profile: layer 12
**2.21 → 1.08**, layer 14 1.43 → 1.16; the total ratio does not move (1.17).

## The beamline (20 Sep)

`simulation/geometry/SND_compact_beamline.xml` is the H2 line as it stood in June
2026: nothing between the vacuum pipe and the box, the box 2–3 m from the vacuum
window — two XCET Cherenkovs (4 Al windows, 0.022 X₀), the exit window (1 mm Al,
0.011 X₀) and 2.5 m of air (0.008 X₀), ≈ 0.04 X₀ in all. Run it with
`COMPACT_FILE=SND_compact_beamline.xml BEAM_Z_MM=-2650`. The HGCAL 2018 budget
of ≈ 0.5 X₀ (JINST 17 P05022), which was tried first, is kept as
`SND_compact_beamline_hgcal.xml` for reference: half of it was HGCAL's own trigger
scintillators, veto and delay wire chambers, and it is not ours. Three knobs
(`Beamline_ResidualAl`, `Beamline_DwcPcb`, `Beamline_ScintThick`, all 1 µm) make a
variant a one-constant edit.

What the scan taught (e⁻ 20 / 52 GeV, th230, per-slab digitizer, chip selection):
layers 0–1 measure the material in front and only that — 1.9 / 1.4 with 0.49 X₀,
1.29 / 1.11 with 0.23, 0.85–0.91 / 0.92–0.94 with the real line; where it sits does
not matter (0.49 X₀ lumped at the box or spread over 1–3 m: 1.87 → 1.75). With the
real line the residual is a ratio rising with depth, 0.85 at layer 0 to 1.3 at
layers 9–14, and a simulated shower 0.3 layer deeper. Split by hit amplitude the
hits ≤ 60 ADC and the hit counts agree in every layer; the excess is entirely in
hits > 300 ADC and grows with depth and with energy — not the trigger, not the
upstream material; candidates are the SKIROC2 response at high occupancy, charge
spreading in the data, and the tungsten thicknesses of the stack (step at layer 9,
where the plates go from 4.2 to 5.6 mm). A Gaussian-core or shape selection of the
data events changes none of this (the data depth distribution is shifted as a
whole, not a tail).

Caveat on the data MIP tables: re-fitting th230 from its merged histograms with
the current calibrator (`calibration/MuonCalib_gaudi_slabwindow/`, not deployed)
reproduces slab 12's 94 masked channels and MPV 39 — the window was never the
limit, the truncation is — and changes ~17 % of the *other* channels' entries
relative to the July table, i.e. the deployed tables are not reproducible with the
current calibrator. Re-deploying means re-reconstructing every data tree; left
for a dedicated pass (a truncated-Landau fit for slab 12 with it).

### The hit selection: `HIT_SELECTION=chip` against `EcalEventBuilder HitSelection=adc`

The data's event builder keeps a channel only if its `hitbit_high` is set
(`bestScaPerChannel`). On the raw chunks, ~15% of the high-gain ADC of every
layer sits in channels of a triggered chip with a clear signal and no hit bit
(80% with no bit in any SCA of the window, 17% only in an earlier SCA, 3% only
in the retrigger SCA), the same at 20 / 52 / 74 GeV and flat over layers 1-12.
`siwecal-tb2026` gained `HitSelection='adc'` (bit in any SCA **or**
`adc_high − pedestal > AdcHitThreshold` = 30 ADC, read at the largest SCA), and
this repo the mirror image:

| env (job3) | RealDigitizer | AdcDigitizer | tree |
|---|---|---|---|
| `HIT_SELECTION=cell` | one hit per cell whose discriminator fired | every hit kept | — |
| `HIT_SELECTION=chip` (default) | every sampled cell is written, with `SiPadHitsRealDigitizedTrigger` (1 = fired); a cell that did not fire is sampled at fast peak + delay, where the chip's hold lands | `TriggerCollection` on: a chip counts when a cell fired; on it a cell is a hit if it fired or `hit_hg > AdcHitThreshold`; verdict in `SiPadHitsRealAdcKept` (still 1:1) | `sim_to_ecal_tree` drops the 0s |

Redone like for like (data trees rebuilt with `adc`, gains rescanned, roll-overs
refitted, Gaussian core μ of the event ADC sum):

| set | point | gain old → new | sim/data, hit bit / cell | sim/data, adc / chip |
|---|---|---|---|---|
| th230 | e⁻ 20 GeV run 20 | 21.9 → 22.2 | 1.16 | **1.14** |
| th230 | e⁻ 52 GeV run 13 | | 1.18 | **1.15** |
| th220 | e⁻ 74 GeV run 72 | 21.5 → 21.6 | 1.36 | **1.27** |
| th210 | e⁻ 52 GeV eudaq 287 | 19.0 → 19.05 | 1.16 | **1.08** |
| th210 | e⁻ 74 GeV eudaq 286 | | 1.28 | **1.13** |

The data recover 15-25%, but so does the simulation where its measured plateau
efficiency is 0.86 (th230, th220): that plateau *was* the missing-bit effect,
applied cell by cell. The gap that remains is a longitudinal one — the sim/data
ratio per layer rises from 0.55 at layer 0 to 1.5 at layer 14 — i.e. the
simulated shower is deeper and its sampling fraction larger: geometry (absorber
thickness/density, upstream material), not electronics. The MIP tables were made
with the hit-bit selection; redo them before using `adc` for calibrated energies.
Results at `Simulation/Processed/adc_vs_tb/final_chip/`.

### Why RealDigitizer runs before GeV2MIPConversion

The shaping needs the per-step `CaloHitContributions`, and
`GeV2MIPConversion` — like `BasicDigitizer`, `DetectorFlipper` and
`ChannelMapper` — creates fresh hits copying only CellID, position and energy.
So the real chain reads `SiPadHits` in GeV and normalises to MIP itself through
`MIPValues`; `GeV2MIPConversion` would be redundant on that branch.

### How the digitised energy and time reach the analysis

`HitEnergyContent='digitized'` (the default) writes the shaped slow-sample
amplitude, in MIP, into the **hit energy**. That is what makes the digitisation
visible downstream: `DetectorFlipper:246`, `ChannelMapper:231` and
`analysis/sim_to_ecal_tree.py` all read `hit.getEnergy()` and nothing else.

The trigger time has nowhere to live in a `SimCalorimeterHit`, so it travels in
the parallel `podio::UserDataCollection<float>` named by
`DigitizedTimeCollection`. That works because **`DetectorFlipper` and
`ChannelMapper` are 1:1** — one output hit per input hit, never dropping any —
so the vector stays index-aligned all the way to `SiPadHitsRealMapped`, the same
way `SiPadHitsMasked` already does. `sim_to_ecal_tree.py` reads it into the
`hit_time` branch, with a length check guarding against pairing hits with the
other chain's times.

```bash
# both chains into one file, then analyse the simple one against the real one
CALIB_THRESHOLD=th230 DIGI_MODE=both bash gaudi_jobs/1_mu_beam_pipeline/1_mu_pipeline.sh
bash analysis/run_pid_sim.sh --collection SiPadHitsMapped --format both
```

### The shaping preserves the amplitude by construction

`crRcResponse` is normalised to **unit peak gain**: one instantaneous step of
A MIP produces a pulse whose maximum is exactly A MIP, at `t = tauNs`, for any
order. That normalisation is `exp(n)/n^n`, and
`analysis/tests/test_cell_shaping.py` pins it by compiling the header.

It used to be `4/n!`, which is neither unit peak nor unit area: at order 2 it sat
**8.3% above** unit peak ((4/2!)·2²·e⁻² = 1.0827), so every digitised energy
carried a ~7% scale on top of the deposit — and the scale changed with the order
(×1.47 at order 1, ×0.90 at order 3) with nothing to announce it.

What is left is deliberate, and measured on the 1000-event samples as
`sum_energy` of the real chain over the simple chain:

| sample | real / simple | hits kept |
|---|---|---|
| mu- 100 GeV | 0.992 | 17.2 of 17.5 |
| e- 52 GeV | 0.986 | 575 of 632 |
| e- 74 GeV | 0.987 | 720 of 793 |

Two effects, both physical: the slow channel is sampled at `triggerTime + DelayNs`
(160 ns) while its peak is at `tauSlowNs` (180 ns), so the sample rides the rise
and lands 1–3% low, by an amount that depends on the amplitude through the
trigger's time walk; and the 0.5 MIP fast-channel threshold drops the smallest
hits. Setting `DelayNs = TauSlowNs` samples the peak and removes the first.

### Where the layer z comes from

There is one z table per frame and nothing may hold a private copy:

| Frame | Source | Values |
|---|---|---|
| Simulation | the compact XML, via `sipad::sensitiveLayers` (`gaudi_source/SiPadLayerGeometry.h`) | 49.35 … 274.275 mm, pitch 15 mm |
| Test beam | `mappings/slab_z_positions.yml` | 0 … −225 mm, pitch 15 mm |

`Ecal_LayerDistance` is 15 mm and layers 10→11 are two pitches apart (the empty
rail slot). The two tables are the same detector, mirrored and offset.

`ACTSGeoSvc` and `DetectorFlipper` both read the simulation frame through the
same helper, so they cannot disagree about where layer N is. `DetectorFlipper`
has **no built-in table**: leaving `ZPositions` unset takes z from DD4hep and
emits a WARNING saying so (that is a no-op flip — the hits keep the simulation
frame). When `ZPositions` *is* given, its layer spacing is cross-checked against
the geometry and a mismatch is warned about.

That guard exists because the pitch went 11 mm → 15 mm in July 2026 and three
separate copies of the table stayed behind — `DetectorFlipper`'s default (on a
third value, 16.6 mm) and two test files. A wrong z table does not fail, it
produces perfectly plausible hits in the wrong place. The tests read the same
two sources through `analysis/tests/geometry_ref.py`.

---

## Job 4 — Tracking

**Files:** `gaudi_source/ShowerTagger.cpp`, `SiPadMeasConverter.cpp`, `ACTSProtoTracker.cpp`, `ACTSGeoSvc.cpp`
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

### Single-file output — no separate tracks.edm4hep.root

The job itself is a plain Gaudi `IOSvc` reader/writer, so it always needs a
distinct `OUTPUT_FILE` name (writing back onto the file it is still reading
would corrupt it). Every calling `.sh` script therefore points `OUTPUT_FILE`
at a temp name and, once `k4run` exits, `mv`s it back onto
`digitized.edm4hep.root`:

```bash
INPUT_FILE="digitized.edm4hep.root" INPUT_COLLECTION="SiPadHitsDigi" \
    OUTPUT_FILE="digitized_tracks_tmp.edm4hep.root" SEED_MOMENTUM=100.0 \
    k4run ../pid2026_common/job4_tracking.py
mv digitized_tracks_tmp.edm4hep.root digitized.edm4hep.root
```

The temp name must still **end in `.root`**: `IOSvc`'s `Writer` picks its
backend off the filename, and a name like `digitized.edm4hep.root.tracks_tmp`
(the first thing tried here) fails initialization with `Unknown file type for
file ... with type default` — caught by actually running the condor chunk
script locally, not just `bash -n`.

`outputCommands = ["keep *"]` means the swapped-in file still carries every
collection `digitized.edm4hep.root` had (`SiPadHitsDigi`, `SiPadHitsFlipped`,
`SiPadHitsMapped`, `SiPadHitsMasked`, …) *plus* `ACTSTracks`, `EMShowers`,
`SiPadShowerFlags` and `SiPadMeasurements` — one edm4hep file staged to
`Processed/` (`<label>_digitized.edm4hep.root`), not two overlapping ones.
This used to write a standalone `tracks.edm4hep.root` (staged separately as
`<label>_tracks.edm4hep.root` in some pipelines, and not staged at all — so
silently dropped — in others); every `1_*_pipeline.sh` and the condor chunk
scripts (`simulation/run_script/generic_condor_beam_chunk.sh`,
`reprocess_chunk.sh`) now use the temp+swap pattern instead. `job5_rntuple.py`
reads `digitized.edm4hep.root`, so `event_display` gets `ACTSTracks` straight
from it (`ShipHits.root`'s `ACTSTracks`/`SiPadMeas` RNTuples).

**The PID chain merges it in too — still one file, no companions.**
`analysis/run_pid_sim.sh` feeds `digitized.edm4hep.root` to
`analysis/sim_to_ecal_tree.py`, which writes `event = <original frame index>`
for every row; `k4SiWEcalReco`'s `EcalToEDM4hep` (in `siwecal-tb2026`) reads
that tree in the same order, so its own frame index is that same original
frame index. `run_pid_sim.sh` then passes that exact `digitized.edm4hep.root`
path **explicitly** to `run_pid_batch.py --tracks-file` (no sibling-file
auto-discovery — the caller already knows the path, so it just says so); `run_
pid_batch.py` (also in `siwecal-tb2026`) forwards it to
`siwecal_common.edm4hep_pid.write_filtered`, which merges `ACTSTracks`/
`EMShowers`/`SiPadMeasurements`/`SiPadShowerFlags` straight into the PID
output at that matching index — so `ecal_sim.edm4hep.root` (or the staged
`<label>_ecal.edm4hep.root`) itself carries the tracks: `event_viewer` reads
them from the one file it already has open
(`PidFileReader.track_counts()` / `Edm4hepEventReader.n_tracks()`), no sibling
lookup anywhere in the chain. `SiPadHits*`/`MCParticles` are deliberately left
out of that merge —
`ECalHits` already covers the raw hits, and those collections'
`SimCalorimeterHit`→`CaloHitContribution` relations segfault when resolved
against a second, concurrently-open podio `Reader`; the four merged
collections carry no such relations. See `siwecal-tb2026/README.md` ("EDM4hep
PID file") and `docs/acts_integration.md` for more.

**Step 1 — `ShowerTagger`:** identifies electromagnetic cascades and keeps
their hits out of the ACTS measurement pool — see below. Writes per-hit veto
flags (`SiPadShowerFlags`) and reconstructed showers (`EMShowers`).

**Step 2 — `SiPadMeasConverter`:** turns the pad hits into
`edm4hep::TrackerHit3D` measurements, writing the layer index into `quality` so
the tracker can look the surface up by address. Variance is `pitch²/12` per
axis. Positions are taken from the hit, not recomputed.

**Step 3 — `ACTSProtoTracker`:** Hough seeding → CKF → KalmanFitter refit →
event-level deduplication, writing `ACTSTracks` with one `AtIP` state carrying
the seed position and one `AtOther` state per surface.

### Showers are not tracks

A track through an electromagnetic cascade is not a physical object: the pads
are secondaries spraying transversely, not samples of a trajectory. Fitting a
line through them yields a well-formed track with a small chi2 and no meaning.
So shower hits are identified **before** the measurement pool is built and never
reach ACTS; what an EM event produces is a *shower*, not a track.

This bites harder here than in a tracker-plus-calorimeter setup: every layer is
1.2-2.0 X0, so a high-energy electron is already showering in layer 0 or 1 and
there is no incoming segment left to fit.

`ShowerTagger` counts hits per layer. A MIP lights one or two pads per layer, a
shower lights tens; `ShowerMinConsecutive` consecutive layers at or above
`ShowerNHitsThreshold` hits mark the onset (two layers, so a single delta-ray
spike does not count). Hits from the onset onwards are flagged. If fewer than
`MinTrackLayers` layers precede the onset, **every** hit is flagged — a stub of
one to three points cannot define a trajectory, and letting it through is how a
shower event ends up with a "track" anyway.

| Property | Default | Description |
|---|---|---|
| `ShowerNHitsThreshold` | 4 | Hits in one layer for that layer to count as dense |
| `ShowerMinConsecutive` | 2 | Consecutive dense layers required to declare an onset |
| `MinTrackLayers` | 4 | Pre-shower layers needed to still offer the segment to the tracker |
| `Enabled` | `True` | `False` disables tagging (all flags zero, no showers) |

The veto reaches ACTS through `SiPadMeasConverter.InputFlags`, the same
per-hit-flag idiom as `ChannelMapper`'s `OutputMaskedFlags`. The flags are
positional, so the converter fails loudly if the two collections have different
lengths rather than vetoing the wrong hits.

`EMShowers` is an `edm4hep::ClusterCollection`: `type = 1`, `energy` = total (in
the input's units, MIPs after digitisation), `position` = energy-weighted
barycentre, and `shapeParameters` = [start layer, layer of maximum, layers
spanned, transverse RMS in mm, number of hits].

Behaviour by particle:

| | Result |
|---|---|
| Muon | No onset; the fifteen-layer track is untouched |
| Electron | Onset at layer 0-1, so a shower and **no track** |
| Radiative muon | Onset late: keeps the incoming track *and* gets a shower |
| Pion punching through then showering | Same — incoming track plus shower, which is what PID wants |

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
| Showers reconstructed | 49 (radiative muons) | 299/300 |
| Events with a track | 987/1000 | **0** |
| Tracks per event | 1.00 | 0 |
| Track spans all 15 layers | 99.8% | — |
| Residual vs hits (rms) | 0.21 mm (x), 0.28 mm (y) | — |
| Hits vetoed as shower | 0% in 951 events | 100% |

Zero tracks in the electron sample is the intended result, not a failure. The 13
muon events without a track are ones where the muon dumped a real EM cascade
early (`ShowerTagger` finds an onset in 49 events; in 12 of them it starts at
layer 2-3, leaving too short a stub). Median reconstructed shower: 4990 MIP,
onset at layer 1, maximum at layer 6, transverse RMS 22.6 mm.

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

## Comparing the digitisation with test-beam data (ADC)

**File:** `analysis/compare_adc_data_sim.py`

The simulation's ecal tree has no ADC at all (`hit_hg` is 0; the digitised
amplitude lives in `hit_energy`, in MIP), while the test-beam tree's `hit_hg`
IS the pedestal-subtracted high-gain ADC. The comparison is made on the data's
own axis by pushing the simulated MIP back through the SAME per-channel MIP
table the data reconstruction used:

```
ADC_sim(hit) = hit_energy[MIP] x mpv(slab, chip, channel)
```

which is only possible because `ChannelMapper` has already rewritten the
simulated CellIDs into the test-beam `slab/chip/channel` format, and masks the
channels that table has no MPV for — so both samples are missing the same ones.

```bash
python3 -m analysis.compare_adc_data_sim \
    --data /eos/.../Reconstruction/TB2026CERN_run_000013/ecal_TB2026CERN_run_000013.root \
    --sim digi=ecal_sim_e52_real.root --sim simple=ecal_sim_e52_simple.root \
    --mip-file ../siwecal-tb2026/calibration/MuonCalib_gaudi/mips/th230/\
MIP_pedestalsubmode1_TB2026CERN_run_000004_highgain.txt \
    --tag e52 --outdir plots/
```

Two data curves are always drawn, because above `AdcSaturationThreshold`
(1500 ADC) the two are not the same measurement:

| curve | what it is |
|---|---|
| `data (raw hit_hg)` | the detector reading, high-gain saturation included |
| `data (linearised)` | `hit_energy x mpv` — the ADC an unsaturated preamp would have given, recovered from the low gain |

The simulation has **no high-gain saturation model**, so `data (linearised)` is
the apples-to-apples curve; the gap between the two data curves above 1500 ADC
is exactly what such a model would have to reproduce.

### What the comparison currently says

Measured on 1000-event samples against the P1 runs at the same energy
(hits above 0.5 MIP, unmasked, calibrated channels only):

| sample | hits/event | <ADC>/hit | event sum [ADC] | sigma/mu | sim/data |
|---|---|---|---|---|---|
| e- 52 GeV data (run 13) | 380 | 149 | 56 600 | 0.30 | 1.00 |
| e- 52 GeV sim (digi) | 569 | 200 | 114 100 | 0.066 | **2.01** |
| e- 52 GeV sim (simple) | 573 | 201 | 115 200 | 0.066 | 2.04 |
| e- 74 GeV data (run 7) | 453 | 147 | 66 600 | 0.37 | 1.00 |
| e- 74 GeV sim (digi) | 712 | 226 | 160 900 | 0.053 | **2.42** |
| e- 74 GeV sim (simple) | 717 | 227 | 162 600 | 0.053 | 2.44 |

Three things to keep in mind before reading that ratio as a digitisation bug:

- **The single-hit ADC scale is right.** On a MIP-like run the simulated Landau
  sits on the data's, peak within ~20% and the tail on top of it, so the
  MIP -> ADC conversion and the per-hit response are not what is off.
- **The discrepancy is in the shower, and grows with energy.** Layers 0-2 agree;
  the data profile then flattens where the simulation peaks. The simulated
  response is linear between the two energies (x1.41 for x1.42 in beam energy),
  the data's is not (x1.18), and the data resolution gets *worse* with energy
  (0.30 -> 0.37) — the signature of something limiting the data at high
  occupancy (SCA depth, retriggers, BCID splitting, event selection), not of the
  shaping.
- **The flat part of the factor is a calibration scale, not saturation.**
  `analysis/plot_response_scan.py` shows the data at 33-37 MIP/GeV from 7.5 to
  52 GeV against the simulation's 69, and a factor that does not move with energy
  cannot be an occupancy effect. The th230 MIP table's MPV is 1.47x th220's,
  channel by channel, while the pedestal means (245 ADC), the pedestal widths
  (1.7-2.4 ADC) and the LG->HG anchor slope (k = 0.0925/0.0963/0.0961) are the
  same across the three threshold sets — so the ADC scale did NOT change with the
  threshold and the MPV difference is a calibration artifact. Only 39.6% of the
  th230 channels carry a real per-channel fit (98.5% at th220); 54% take the
  chip-level fallback (`empv = -3`), on a median of 183 entries per channel.
  `analysis/mip_threshold_bias.py` measures what that costs: pushing the
  *simulated* muon sample through a threshold and recovering the MPV the way the
  calibration does, a bias of x1.47 needs a threshold at 1.35x the MIP MPV, which
  keeps 40% of the MIP hits.
- **The data is unselected.** Both runs are read event by event with no beam or
  quality cut, as the validation does (`kept_frac = 1`).

---

## Fast shaper vs slow shaper

**File:** `analysis/compare_shapers.py`

The two sides record the fast channel differently, so they cannot be overlaid
directly:

| | fast shaper (trigger) | slow shaper (amplitude) |
|---|---|---|
| simulation | peak [MIP] -> `hit_fast`, compared with `Threshold` | sample at `triggerTime + DelayNs` -> `hit_energy` |
| test beam | `hitbit_high`, one BIT per channel and SCA | `adc_high` -> `hit_energy` |

What is comparable is the **turn-on**: the probability that a cell enters the
event at all, against the amplitude the slow shaper measured for it. On the data
side that is `P(hitbit_high = 1 | amplitude)` — and it *is* the data's hit
selection, because `EventBuilder::bestScaPerChannel` drops every channel whose
bit never fired. On the simulation side it is obtained by matching the `simple`
chain (Threshold = 0, so it keeps every cell) against the `real` one, channel by
channel.

The fast peak reaches the analysis through `SiPadHitsRealDigitizedFast`, a
parallel `UserDataCollection<float>` written by `RealDigitizer` next to the
digitised energy and time, which `sim_to_ecal_tree.py` puts in the `hit_fast`
branch.

### Measured on run 13 (th230) against the 52 GeV sample

| | data | sim |
|---|---|---|
| turn-on 50% | 0.80 MIP | 0.53 MIP |
| turn-on 90% | — (never reached) | 0.63 MIP |
| plateau efficiency above 2.5 MIP | **0.868** | 1.000 |

**The simulation's trigger is a step at 0.5 MIP that is 100% efficient; the real
one turns on at 0.80 MIP, takes ~0.5 MIP to get there, and plateaus at 87%.**
`Threshold` is also a single hardcoded number, so nothing about it changes when
the sample is meant to match th210 or th220 — the threshold is the one thing in
the chain that genuinely is per-threshold-set.

What that costs, replaying the measured turn-on over the simulated hits
(e- 52 GeV, per event):

| scenario | hits | MIP | sim/data |
|---|---|---|---|
| sim, current step at 0.5 MIP | 569 | 3579 | 2.01 |
| sim + measured turn-on | 472 | 3093 | 1.74 |
| sim + measured turn-on + MPV corrected x1.47 | 425 | 3043 | **1.16** |

(data: 380 hits, 1780 MIP as calibrated today, 2617 MIP if the th230 MPV is
indeed x1.47 high. The two corrections are not independent — a biased MPV also
moves the measured turn-on along its own amplitude axis — so they are applied
together, not multiplied.)

Two independently motivated fixes take the discrepancy from x2.0 to x1.16. The
shaping itself is not implicated either way: `slow / fast` per hit has a median
of 0.991 with a 3.5% IQR, and the residual trend with amplitude (1.047 at
0.5 MIP to 0.988 at 80 MIP) is the time walk of sampling at
`triggerTime + 160 ns` when the slow peak is at 180 ns.

---

## Adding a New Gaudi Algorithm

1. Create `gaudi_source/MyAlgorithm.cpp` with `DECLARE_COMPONENT(MyAlgorithm)` at end
2. Add `.cpp` to `gaudi_add_module(SND_reco SOURCES ...)` in `CMakeLists.txt`
3. `/build` to rebuild
4. Import in Python: `from Configurables import MyAlgorithm`

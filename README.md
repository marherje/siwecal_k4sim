# siwecal_k4sim

Full key4hep simulation + reconstruction chain for the SiW-ECAL test beam 2026.

**Stack:** DD4hep v01-35 / key4hep (release pinned in `.key4hep-release`) / Gaudi v40 / ROOT 6.38 / AlmaLinux 9 (lxplus)

**Companion repo:** `../siwecal-tb2026` — k4SiWEcalReco plugin, calibration

---

## Quick start (daily use)

```bash
cd /path/to/siwecal_k4sim
source init_siwecal_soft.sh   # loads everything: key4hep + local build + event viewer
```

After sourcing, the following are all available in the same shell:

| Command | What it does |
|---------|-------------|
| `k4run <job>.py` | Run Gaudi digitisation / reconstruction jobs |
| `ddsim ...` | Run Geant4 simulation via DD4hep |
| `bash analysis/run_pid_sim.sh` | Convert sim → ecal TTree → k4SiWEcalReco shower vars |
| `python event_display/event_display_eve.py` | 3D TEve event display |
| `bash event_viewer/launch_sim.sh` | Dash web event viewer (http://localhost:8050) |

---

## First-time setup (run once per account/machine)

### Prerequisites

- Access to CVMFS (`/cvmfs/sw.hsf.org/key4hep/`)
- `siwecal-tb2026` checked out **next to** this repo:
  ```
  parent/
  ├── siwecal_k4sim/     ← this repo
  └── siwecal-tb2026/    ← companion repo
  ```

### Step 1 — Build the Gaudi algorithms

```bash
source init_siwecal_soft.sh   # loads key4hep
bash build.sh                  # cmake + make install → install/
```

> `build.sh` runs cmake with `RelWithDebInfo` and installs into `install/`.
> Repeat this step any time a `.cpp` file in `gaudi_source/` changes.

### Step 2 — Create the event_viewer virtualenv & check k4SiWEcalReco

```bash
bash setup_venv.sh
```

This script creates `.venv-viewer` (in this repo root) with `dash` and `plotly`
on top of key4hep Python.

> After this step, `source init_siwecal_soft.sh` is all you need every session.

---

## Repository layout

```
siwecal_k4sim/
├── init_siwecal_soft.sh      ← daily env setup (source this)
├── setup_venv.sh             ← one-time setup (run once)
├── build.sh                  ← build Gaudi algorithms
├── CMakeLists.txt
│
├── gaudi_source/             ← C++ Gaudi algorithms
│   ├── GeV2MIPConversion.cpp
│   ├── BasicDigitizer.cpp
│   ├── DetectorFlipper.cpp   ← rewrite z positions per layer
│   ├── ChannelMapper.cpp     ← remap CellID sim→TB + MIP masking
│   └── ...
│
├── gaudi_jobs/               ← per-pipeline job scripts
│   ├── 1_mu_beam_pipeline/   ← muon beam (5 GeV)
│   ├── 1_e_beam_pipeline/    ← electron beam (5 GeV)
│   ├── 1_e54_beam_pipeline/  ← electron beam (54 GeV, centre −45/+45 mm)
│   ├── 1_mu_PG_pipeline/     ← muon particle gun (54 GeV)
│   ├── 1_e_PG_pipeline/      ← electron particle gun (54 GeV)
│   └── mip_extraction_pipeline/
│
├── simulation/
│   ├── geometry/             ← DD4hep compact XML
│   └── run_script/           ← HTCondor submission scripts
│       ├── launch_PG.sh      ← submit particle-gun jobs
│       ├── launch_beam.sh    ← submit Gaussian-beam jobs
│       ├── generic_condor_PG.sh
│       └── generic_condor_beam.sh
│
├── analysis/
│   ├── sim_to_ecal_tree.py   ← digitized.edm4hep.root → ecal TTree (hit_X0, ...)
│   ├── run_pid_sim.sh        ← sim→ecal tree + k4SiWEcalReco in one step
│   └── tests/
│
├── mappings/                 ← pad maps, slab-z positions, W thicknesses
│
├── masking_info/
│   └── calibration/          ← MuonCalib_gaudi/mips/th{210,220,230}/ — the MIP
│                               tables ChannelMapper masks on. dummy_mip_map_*
│                               masks nothing, for A/B checks. Pedestals, LG→HG
│                               anchors and .diagnostics.root are not tracked:
│                               nothing here reads them, see ../siwecal_calib_archive
│
├── event_display/            ← 3D TEve display (ROOT-based, works offline)
│   ├── event_display_eve.py
│   ├── detector_config.json
│   └── launch.sh
│
└── event_viewer/             ← Dash web viewer (self-contained)
    ├── launch_sim.sh
    ├── sim_settings.yml
    └── [Python package — app.py, model/, ui/, viz/, ...]
```

---

## Running a pipeline

All pipelines follow the same pattern. Example for the 54 GeV electron beam:

```bash
source init_siwecal_soft.sh
bash gaudi_jobs/1_e54_beam_pipeline/1_e54_beam_pipeline.sh
```

The pipeline steps (simplified — no shuffler/splitter/tracker needed for direct sim):

```
simulation output (.edm4hep.root)
  → [job3: GeV2MIP + BasicDigitizer + DetectorFlipper + ChannelMapper]
  → digitized.edm4hep.root   (SiPadHitsMapped, SiPadHitsMasked)
  → [sim_to_ecal_tree.py]
  → ecal_sim.root             (TTree: hit_z, hit_X0, chip, chan, energy, ...)
  → [k4SiWEcalReco]
  → ecal_sim.edm4hep.root     (ECalHits, ECalPid with barz/moliere/...)
  → ecal_sim.valtree.root     (flat TTree for analysis)
```

### ecal TTree variables

| Branch | Type | Description |
|--------|------|-------------|
| `hit_slab` | `int[nhit_chan]` | Layer index (0–14) |
| `hit_chip` | `int[nhit_chan]` | Chip ID (0–15) |
| `hit_chan` | `int[nhit_chan]` | Channel ID (0–63) |
| `hit_energy` | `float[nhit_chan]` | Hit energy [MIP] |
| `hit_x/y/z` | `float[nhit_chan]` | Hit position [mm] |
| `hit_X0` | `float[nhit_chan]` | Cumulative W radiation length [X₀] |
| `hit_ismasked` | `int[nhit_chan]` | 1 = masked/uncalibrated channel |

`hit_X0` values: 0.80 (layer 0, 2.8 mm W) → 20.40 X₀ (layer 14, 71.4 mm W total).

---

## Launching the event display (TEve / ROOT)

Reads `ShipHits.root` from the RNTuple pipeline (job5_rntuple.py).

```bash
source init_siwecal_soft.sh
cd event_display
python event_display_eve.py --hits ../gaudi_jobs/1_mu_beam_pipeline/ShipHits.root --window 0
```

---

## Launching the event viewer (Dash web app)

Reads `ecal_sim.edm4hep.root` or `ecal_sim.valtree.root` produced by `run_pid_sim.sh`.

```bash
source init_siwecal_soft.sh
bash event_viewer/launch_sim.sh                                          # mu-beam default
bash event_viewer/launch_sim.sh --data-dir gaudi_jobs/1_e54_beam_pipeline  # e- 54 GeV
bash event_viewer/launch_sim.sh --port 8051 --debug
```

Then open **http://localhost:8050** in a browser. On lxplus, forward the port with SSH:

```bash
ssh -L 8050:localhost:8050 lxplus.cern.ch
```

---

## Submitting simulations on HTCondor

```bash
source init_siwecal_soft.sh
cd simulation/run_script

# Particle-gun (monoenergetic pencil beam)
bash launch_PG.sh          # mu- and e- at 5, 50, 54 GeV, several positions

# Gaussian beam (realistic beam profile)
bash launch_beam.sh        # mu- and e- at 54 GeV, centre (−45,+45) mm, σ=(20.5,16.5) mm
```

Monitor jobs: `condor_q $USER`

Output files land in `simulation/run_script/data/`.

### PID campaign 2026 (e⁻/μ⁻/π⁻, 74 & 99 GeV, chunked)

50 000 events per (particle, energy) at beam centre **(1.139, 1.164) mm**, split
into 25 condor jobs of 2 000 events each (a single 50k electron job would be
~26 h). Every chunk job runs **ddsim → digitisation → ecal-tree conversion** on
the worker's scratch and stages out only:

```
Generated/chunks/output_<point>_cNNN.edm4hep.root   raw sim (kept for reprocessing)
Processed/chunks/<point>_cNNN_ecal.root             ecal TTree (run = E*1000 + chunk)
```

```bash
cd simulation/run_script
bash launch_beam_pid2026.sh                                   # all 150 jobs
DRY_RUN=1 bash launch_beam_pid2026.sh                         # steer files only
PARTICLES="pi-" ENERGIES="99" CHUNKS="3 17" bash launch_beam_pid2026.sh  # resubmit
```

Beam optics per species: e⁻ σ=(13.75, 8.25), μ⁻ σ=(38.5, 46.75),
π⁻ σ=(24.75, 13.75) mm; σ_E = 2 % for all. Each chunk gets its own RNG seed
(`species offset + E*1000 + chunk`), so the 25 chunks are statistically independent.

Once the chunks are on EOS, nine pipelines turn them into physics samples —
one per particle for 74 GeV, 99 GeV and the merged 74+99 GeV sample:

```bash
bash gaudi_jobs/2_e74_pid_pipeline/2_e74_pid_pipeline.sh
bash gaudi_jobs/2_mu99_pid_pipeline/2_mu99_pid_pipeline.sh
bash gaudi_jobs/2_pimerged_pid_pipeline/2_pimerged_pid_pipeline.sh
# ... {e,mu,pi} x {74,99,merged}; add --allow-partial to run with missing chunks
```

They all delegate to `gaudi_jobs/pid2026_common/pipeline_common.sh`, which
hadds the chunk trees (for `merged`: the chunks of *both* energies — merging
happens at the ecal-tree level, no reprocessing), runs k4SiWEcalReco and writes
`Processed/<label>_ecal.{root,edm4hep.root,valtree.root}`. Work goes through
`/tmp/$USER/siwecal_pid2026`, never AFS — the merged trees are ~1 GB.

#### Chunk housekeeping in Processed/

`Processed/chunks/` is a staging area, not a deliverable: once a particle's
three samples are in `Processed/`, its chunk trees are dead weight. The drivers
handle that — no pipeline deletes chunks on its own, because the 74 GeV chunks
are consumed twice (by the 74 GeV pipeline *and* by the merged one):

```bash
bash gaudi_jobs/pid2026_common/run_particle.sh --particle e-   # 74, 99, merged, then clean
bash gaudi_jobs/pid2026_common/run_all.sh                      # all three particles
bash gaudi_jobs/pid2026_common/run_all.sh --keep-chunks        # process, keep chunks
```

`run_particle.sh` removes that particle's chunk trees from `Processed/chunks/`
only after verifying all three final samples exist and are non-empty; if any
pipeline fails the chunks stay for the retry. The raw simulation chunks in
`Generated/chunks/` are **never** deleted — they are the reprocessable source,
worth ~26 h of CPU per point.

---

## Running tests

```bash
source init_siwecal_soft.sh
pytest analysis/tests/ -v
```

105 tests covering ChannelMapper, ecal tree converter, and geometry mapping.

# siwecal_k4sim

Full key4hep simulation + reconstruction chain for the SiW-ECAL test beam 2026.

**Stack:** DD4hep v01-35 / key4hep 2026-02-01 / Gaudi v40 / ROOT 6.38 / AlmaLinux 9 (lxplus)

**Companion repo:** `../siwecal-tb2026` — k4SiWEcalReco plugin, event_viewer, calibration

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

This script:
1. Creates `../siwecal-tb2026/.venv-viewer` with `dash` and `plotly` on top of key4hep Python
2. Builds `../siwecal-tb2026/k4SiWEcalReco` if not already built

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
├── masking_info/
│   ├── geometry/             ← FEV10/FEV11 chip-channel→(x,y) pad maps
│   └── calibration/          ← MIP calibration files (dummy + MuonCalib_it2)
│
├── event_display/            ← 3D TEve display (ROOT-based, works offline)
│   ├── event_display_eve.py
│   ├── detector_config.json
│   └── launch.sh
│
└── event_viewer/             ← Dash web viewer (uses siwecal-tb2026)
    ├── launch_sim.sh
    └── sim_settings.yml
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

---

## Running tests

```bash
source init_siwecal_soft.sh
pytest analysis/tests/ -v
```

105 tests covering ChannelMapper, ecal tree converter, and geometry mapping.

#!/bin/bash
# Pipeline for 54 GeV electron beam simulation
# Beam centre: (-45, 45) mm  |  sigma: (20.5, 16.5) mm  |  sigma_E: 2%
# Skips shuffler/splitter/tracker: reads sim output directly.
set -uo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
cd "${REPO_ROOT}/gaudi_jobs/1_e54_beam_pipeline"

SIM_FILE="/eos/experiment/drdcalo/siw-ecal/TB2026-06/Simulation/output_beam_e-_54GeV_xy_-45_45_sigx20.5_sigy16.5_sigE0.02.edm4hep.root"

echo "=== Step 1: digitize + flip + channel mapping ==="
INPUT_FILE="${SIM_FILE}" k4run job3_digitize.py 2>&1 | grep -v "^TCling::LoadPCM"

echo "=== Step 2: ecal tree + shower variables ==="
cd "${REPO_ROOT}"
bash analysis/run_pid_sim.sh \
    --input   "${REPO_ROOT}/gaudi_jobs/1_e54_beam_pipeline/digitized.edm4hep.root" \
    --ecal-tree "${REPO_ROOT}/gaudi_jobs/1_e54_beam_pipeline/ecal_sim.root" \
    --outdir  "${REPO_ROOT}/gaudi_jobs/1_e54_beam_pipeline" \
    --format  both 2>&1 | grep -v "^TCling::LoadPCM"

#!/bin/bash
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
cd "${REPO_ROOT}/gaudi_jobs/1_mu_PG_pipeline"

echo "=== Step 1: shuffle ==="
k4run job1_shuffler.py

echo "=== Step 2: split into time windows ==="
k4run job2_splitter.py

echo "=== Step 3: digitize + flip + channel mapping ==="
INPUT_FILE=timewindows.edm4hep.root k4run job3_digitize.py

echo "=== Step 4: tracking ==="
k4run job4_tracking.py

echo "=== Step 5: RNTuple conversion ==="
k4run job5_rntuple.py

echo "=== Step 6: ecal tree + shower variables ==="
cd "${REPO_ROOT}"
bash analysis/run_pid_sim.sh \
    --input   "${REPO_ROOT}/gaudi_jobs/1_mu_PG_pipeline/digitized.edm4hep.root" \
    --ecal-tree "${REPO_ROOT}/gaudi_jobs/1_mu_PG_pipeline/ecal_sim.root" \
    --outdir  "${REPO_ROOT}/gaudi_jobs/1_mu_PG_pipeline" \
    --format  both

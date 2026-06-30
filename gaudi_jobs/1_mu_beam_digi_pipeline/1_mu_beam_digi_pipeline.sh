#!/bin/bash
set -euo pipefail

EOS_BASE="/eos/experiment/drdcalo/siw-ecal/TB2026-06/Simulation"
SIM_FILE="${EOS_BASE}/Generated/output_beam_mu-_50GeV_xy_1_1_sigx23.5_sigy29.7_sigE0.02.edm4hep.root"
LABEL=$(basename "${SIM_FILE}" .edm4hep.root | sed 's/^output_//')
PROCESSED="${EOS_BASE}/Processed"

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
cd "${REPO_ROOT}/gaudi_jobs/1_mu_beam_digi_pipeline"

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

PIPELINE_DIR="${REPO_ROOT}/gaudi_jobs/1_mu_beam_digi_pipeline"
echo "=== Copying outputs to Processed ==="
cp "${PIPELINE_DIR}/digitized.edm4hep.root" "${PROCESSED}/${LABEL}_digitized.edm4hep.root"
[[ -f "${PIPELINE_DIR}/ShipHits.root" ]] && \
    cp "${PIPELINE_DIR}/ShipHits.root" "${PROCESSED}/${LABEL}_ShipHits.root"
echo "=== Done. Outputs in ${PROCESSED}/ ==="
ls -lh "${PROCESSED}/${LABEL}"* 2>/dev/null || true

#!/bin/bash
set -euo pipefail

EOS_BASE="/eos/experiment/drdcalo/siw-ecal/TB2026-06/Simulation"
SIM_FILE="${EOS_BASE}/Generated/output_mu-_E5.edm4hep.root"
LABEL=$(basename "${SIM_FILE}" .edm4hep.root | sed 's/^output_//')
PROCESSED="${EOS_BASE}/Processed"

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
cd "${REPO_ROOT}/gaudi_jobs/1_mu_PG_digi_pipeline"

echo "=== Step 1: shuffle ==="
k4run job1_shuffler.py

echo "=== Step 2: split into time windows ==="
k4run job2_splitter.py

echo "=== Step 3: digitize + flip + channel mapping ==="
INPUT_FILE=timewindows.edm4hep.root k4run job3_digitize.py

echo "=== Step 4: tracking ==="
# Tracks on the PRE-flip collection from digitized.edm4hep.root (job3's
# output), and the result is written back into that same file (temp + swap)
# rather than a separate tracks.edm4hep.root -- ACTSTracks/EMShowers/
# SiPadMeasurements end up in the one edm4hep file that gets staged.
INPUT_FILE="digitized.edm4hep.root" INPUT_COLLECTION="SiPadHitsWindowed" \
    OUTPUT_FILE="digitized.edm4hep.root.tracks_tmp" \
    k4run ../pid2026_common/job4_tracking.py
mv digitized.edm4hep.root.tracks_tmp digitized.edm4hep.root

echo "=== Step 5: RNTuple conversion ==="
k4run job5_rntuple.py

PIPELINE_DIR="${REPO_ROOT}/gaudi_jobs/1_mu_PG_digi_pipeline"
echo "=== Moving outputs to Processed (EOS only) ==="
mv "${PIPELINE_DIR}/digitized.edm4hep.root" "${PROCESSED}/${LABEL}_digitized.edm4hep.root"
[[ -f "${PIPELINE_DIR}/ShipHits.root" ]] && \
    mv "${PIPELINE_DIR}/ShipHits.root" "${PROCESSED}/${LABEL}_ShipHits.root"
echo "=== Done. Outputs in ${PROCESSED}/ ==="
ls -lh "${PROCESSED}/${LABEL}"* 2>/dev/null || true

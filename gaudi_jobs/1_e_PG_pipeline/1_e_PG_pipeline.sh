#!/bin/bash
set -euo pipefail

EOS_BASE="/eos/experiment/drdcalo/siw-ecal/TB2026-06/Simulation"
SIM_FILE="${EOS_BASE}/Generated/output_PG_e-_xyz_1_1_-1000_dir_0_0_1_E54.edm4hep.root"
LABEL=$(basename "${SIM_FILE}" .edm4hep.root | sed 's/^output_//')
PROCESSED="${EOS_BASE}/Processed"

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
cd "${REPO_ROOT}/gaudi_jobs/1_e_PG_pipeline"

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
    OUTPUT_FILE="digitized_tracks_tmp.edm4hep.root" \
    k4run ../pid2026_common/job4_tracking.py
mv digitized_tracks_tmp.edm4hep.root digitized.edm4hep.root

echo "=== Step 5: RNTuple conversion ==="
k4run job5_rntuple.py

echo "=== Step 6: ecal tree + shower variables ==="
cd "${REPO_ROOT}"
bash analysis/run_pid_sim.sh \
    --input   "${REPO_ROOT}/gaudi_jobs/1_e_PG_pipeline/digitized.edm4hep.root" \
    --ecal-tree "${REPO_ROOT}/gaudi_jobs/1_e_PG_pipeline/ecal_sim.root" \
    --outdir  "${REPO_ROOT}/gaudi_jobs/1_e_PG_pipeline" \
    --format  both

PIPELINE_DIR="${REPO_ROOT}/gaudi_jobs/1_e_PG_pipeline"
echo "=== Moving outputs to Processed (EOS only) ==="
mv "${PIPELINE_DIR}/digitized.edm4hep.root" "${PROCESSED}/${LABEL}_digitized.edm4hep.root"
for ext in root edm4hep.root valtree.root; do
    src="${PIPELINE_DIR}/ecal_sim.${ext}"
    [[ -f "${src}" ]] && mv "${src}" "${PROCESSED}/${LABEL}_ecal.${ext}"
done
[[ -f "${PIPELINE_DIR}/ShipHits.root" ]] && \
    mv "${PIPELINE_DIR}/ShipHits.root" "${PROCESSED}/${LABEL}_ShipHits.root"
echo "=== Done. Outputs in ${PROCESSED}/ ==="
ls -lh "${PROCESSED}/${LABEL}"* 2>/dev/null || true

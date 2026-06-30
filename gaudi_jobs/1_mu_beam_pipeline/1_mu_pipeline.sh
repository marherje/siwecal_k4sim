#!/bin/bash
# siwecal_k4sim muon beam pipeline.
# Reads from Generated/, writes final outputs to Processed/.
#
# Run from the repo root:
#   source init_key4hep.sh
#   bash gaudi_jobs/1_mu_beam_pipeline/1_mu_pipeline.sh

EOS_BASE="/eos/experiment/drdcalo/siw-ecal/TB2026-06/Simulation"
SIM_FILE="${EOS_BASE}/Generated/output_beam_mu-_100GeV_xy_1_1_sigx24.75_sigy13.75_sigE0.02.edm4hep.root"
LABEL=$(basename "${SIM_FILE}" .edm4hep.root | sed 's/^output_//')
PROCESSED="${EOS_BASE}/Processed"

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
cd "${REPO_ROOT}/gaudi_jobs/1_mu_beam_pipeline"

echo "=== Input:  ${SIM_FILE} ==="
echo "=== Label:  ${LABEL} ==="
echo "=== Output: ${PROCESSED}/ ==="

echo "=== Step 1: digitize + flip + channel mapping ==="
INPUT_FILE="${SIM_FILE}" k4run job3_digitize.py

echo "=== Step 2: ecal tree + shower variables ==="
cd "${REPO_ROOT}"
bash analysis/run_pid_sim.sh --format both

echo "=== Copying outputs to Processed ==="
cp "${REPO_ROOT}/gaudi_jobs/1_mu_beam_pipeline/digitized.edm4hep.root" \
   "${PROCESSED}/${LABEL}_digitized.edm4hep.root"
for ext in edm4hep.root valtree.root; do
    src="${REPO_ROOT}/gaudi_jobs/1_mu_beam_pipeline/ecal_sim.${ext}"
    [[ -f "${src}" ]] && cp "${src}" "${PROCESSED}/${LABEL}_ecal.${ext}"
done

echo "=== Done. Outputs in ${PROCESSED}/ ==="
ls -lh "${PROCESSED}/${LABEL}"* 2>/dev/null || true

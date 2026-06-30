#!/bin/bash
# Pipeline for 54 GeV electron beam simulation
# Beam centre: (-45, 45) mm  |  sigma: (13.75, 8.25) mm  |  sigma_E: 2%
# Reads from Generated/, writes final outputs to Processed/.
set -uo pipefail

EOS_BASE="/eos/experiment/drdcalo/siw-ecal/TB2026-06/Simulation"
SIM_FILE="${EOS_BASE}/Generated/output_beam_e-_54GeV_xy_-45_45_sigx13.75_sigy8.25_sigE0.02.edm4hep.root"
LABEL=$(basename "${SIM_FILE}" .edm4hep.root | sed 's/^output_//')
PROCESSED="${EOS_BASE}/Processed"

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
cd "${REPO_ROOT}/gaudi_jobs/1_e54_beam_pipeline"

echo "=== Input:  ${SIM_FILE} ==="
echo "=== Label:  ${LABEL} ==="
echo "=== Output: ${PROCESSED}/ ==="

echo "=== Step 1: digitize + flip + channel mapping ==="
INPUT_FILE="${SIM_FILE}" k4run job3_digitize.py 2>&1 | grep -v "^TCling::LoadPCM"

echo "=== Step 2: ecal tree + shower variables ==="
cd "${REPO_ROOT}"
bash analysis/run_pid_sim.sh \
    --input     "${REPO_ROOT}/gaudi_jobs/1_e54_beam_pipeline/digitized.edm4hep.root" \
    --ecal-tree "${REPO_ROOT}/gaudi_jobs/1_e54_beam_pipeline/ecal_sim.root" \
    --outdir    "${REPO_ROOT}/gaudi_jobs/1_e54_beam_pipeline" \
    --format    both 2>&1 | grep -v "^TCling::LoadPCM"

echo "=== Copying outputs to Processed ==="
cp "${REPO_ROOT}/gaudi_jobs/1_e54_beam_pipeline/digitized.edm4hep.root" \
   "${PROCESSED}/${LABEL}_digitized.edm4hep.root"
for ext in edm4hep.root valtree.root; do
    src="${REPO_ROOT}/gaudi_jobs/1_e54_beam_pipeline/ecal_sim.${ext}"
    [[ -f "${src}" ]] && cp "${src}" "${PROCESSED}/${LABEL}_ecal.${ext}"
done

echo "=== Done. Outputs in ${PROCESSED}/ ==="
ls -lh "${PROCESSED}/${LABEL}"* 2>/dev/null || true

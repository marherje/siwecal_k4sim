#!/bin/bash
# siwecal_k4sim muon beam pipeline.
# Reads from Generated/, writes final outputs to Processed/.
#
# Run from the repo root:
#   source init_key4hep.sh
#   bash gaudi_jobs/1_mu_beam_pipeline/1_mu_pipeline.sh

EOS_BASE="/eos/experiment/drdcalo/siw-ecal/TB2026-06/Simulation"
# The sample in Generated/; the previous name (sigx24.75_sigy13.75) no longer
# exists there. This beam is wider than the 180 mm acceptance, so ~7% of the
# muons miss the detector entirely and come out as empty events -- expected, and
# the reason several tests skip empty events rather than failing on them.
SIM_FILE="${EOS_BASE}/Generated/output_beam_mu-_100GeV_xy_1_1_sigx38.5_sigy46.75_sigE0.02.edm4hep.root"
LABEL=$(basename "${SIM_FILE}" .edm4hep.root | sed 's/^output_//')
PROCESSED="${EOS_BASE}/Processed"

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
cd "${REPO_ROOT}/gaudi_jobs/1_mu_beam_pipeline"

echo "=== Input:  ${SIM_FILE} ==="
echo "=== Label:  ${LABEL} ==="
echo "=== Output: ${PROCESSED}/ ==="

echo "=== Step 1: digitize + flip + channel mapping ==="
INPUT_FILE="${SIM_FILE}" k4run job3_digitize.py

echo "=== Step 1b: ACTS tracking ==="
# Tracks on the PRE-flip collection: DetectorFlipper moves the hit z into
# the test-beam frame, which no longer matches the ACTS surfaces.
# Written to a temp file and swapped back onto digitized.edm4hep.root itself
# (keep * carries the digitised collections forward) so ACTSTracks/EMShowers/
# SiPadMeasurements end up in the ONE edm4hep file that gets staged, instead of
# a separate tracks.edm4hep.root product.
INPUT_FILE="digitized.edm4hep.root" INPUT_COLLECTION="SiPadHitsDigi" \
    OUTPUT_FILE="digitized.edm4hep.root.tracks_tmp" SEED_MOMENTUM=100.0 \
    k4run ../pid2026_common/job4_tracking.py 2>&1 | grep -v "^TCling::LoadPCM"
mv digitized.edm4hep.root.tracks_tmp digitized.edm4hep.root

echo "=== Step 2: ecal tree + shower variables ==="
cd "${REPO_ROOT}"
bash analysis/run_pid_sim.sh --format both

echo "=== Moving outputs to Processed (EOS only) ==="
mv "${REPO_ROOT}/gaudi_jobs/1_mu_beam_pipeline/digitized.edm4hep.root" \
   "${PROCESSED}/${LABEL}_digitized.edm4hep.root"
for ext in root edm4hep.root valtree.root; do
    src="${REPO_ROOT}/gaudi_jobs/1_mu_beam_pipeline/ecal_sim.${ext}"
    [[ -f "${src}" ]] && mv "${src}" "${PROCESSED}/${LABEL}_ecal.${ext}"
done

echo "=== Done. Outputs in ${PROCESSED}/ ==="
ls -lh "${PROCESSED}/${LABEL}"* 2>/dev/null || true

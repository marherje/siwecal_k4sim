#!/bin/bash
# ---------------------------------------------------------------------------
# run_pid_sim.sh
#
# Step 6 of the siwecal_k4sim pipeline:
#   1. Convert digitized.edm4hep.root  →  ecal_sim.root  (ecal TTree format)
#   2. Run k4SiWEcalReco (from siwecal-tb2026) on ecal_sim.root
#      to compute shower variables and write ecal_sim.edm4hep.root /
#      ecal_sim.valtree.root for the event viewer and validation.
#
# Prerequisites
# -------------
#   a) key4hep environment sourced  (source init_key4hep.sh)
#   b) k4SiWEcalReco built in siwecal-tb2026:
#        cd ../siwecal-tb2026
#        source setup.sh
#        cmake -S k4SiWEcalReco -B k4SiWEcalReco/build && \
#          cmake --build k4SiWEcalReco/build -j$(nproc)
#
# Usage
# -----
#   cd /path/to/siwecal_k4sim
#   source init_key4hep.sh
#   bash analysis/run_pid_sim.sh [--validation] [--format {edm4hep,valtree,both}]
#                                [--max-events N] [--outdir DIR]
#
# Outputs (default: gaudi_jobs/1_mu_beam_pipeline/)
#   ecal_sim.root           ecal TTree (intermediate, kept for inspection)
#   ecal_sim.edm4hep.root   k4SiWEcalReco output with shower variables
#   ecal_sim.valtree.root   plain TTree variant (with --format valtree/both)
# ---------------------------------------------------------------------------

set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
TB2026_ROOT="$(cd "${REPO_ROOT}/../siwecal-tb2026" 2>/dev/null && pwd)" || {
    echo "ERROR: siwecal-tb2026 not found at ${REPO_ROOT}/../siwecal-tb2026"
    echo "  Clone it next to siwecal_k4sim and try again."
    exit 1
}

# --------------------------------------------------------------------------- #
# Defaults (can be overridden via CLI flags below)
# --------------------------------------------------------------------------- #
INPUT_EDM4HEP="${REPO_ROOT}/gaudi_jobs/1_mu_beam_pipeline/digitized.edm4hep.root"
ECAL_TREE="${REPO_ROOT}/gaudi_jobs/1_mu_beam_pipeline/ecal_sim.root"
OUT_DIR="${REPO_ROOT}/gaudi_jobs/1_mu_beam_pipeline"
COLLECTION="SiPadHitsMapped"
FORMAT="edm4hep"
VALIDATION_FLAG=""
MAX_EVENTS=""

# --------------------------------------------------------------------------- #
# Parse arguments
# --------------------------------------------------------------------------- #
while [[ $# -gt 0 ]]; do
    case "$1" in
        --validation)       VALIDATION_FLAG="--validation"; shift ;;
        --format)           FORMAT="$2"; shift 2 ;;
        --format=*)         FORMAT="${1#*=}"; shift ;;
        --max-events)       MAX_EVENTS="--max-events $2"; shift 2 ;;
        --max-events=*)     MAX_EVENTS="--max-events ${1#*=}"; shift ;;
        --outdir)           OUT_DIR="$2"; shift 2 ;;
        --outdir=*)         OUT_DIR="${1#*=}"; shift ;;
        --input)            INPUT_EDM4HEP="$2"; shift 2 ;;
        --input=*)          INPUT_EDM4HEP="${1#*=}"; shift ;;
        --ecal-tree)        ECAL_TREE="$2"; shift 2 ;;
        --ecal-tree=*)      ECAL_TREE="${1#*=}"; shift ;;
        -h|--help)
            sed -n '/^# Usage/,/^# Outputs/p' "$0" | grep -v '^#---'
            exit 0 ;;
        *)
            echo "Unknown flag: $1  (use --help for usage)"; exit 1 ;;
    esac
done

# --------------------------------------------------------------------------- #
# Verify k4SiWEcalReco build
# --------------------------------------------------------------------------- #
K4RECO_BUILD="${TB2026_ROOT}/k4SiWEcalReco/build"
if [[ ! -f "${K4RECO_BUILD}/libk4SiWEcalRecoPlugins.so" ]]; then
    echo "ERROR: k4SiWEcalReco not built. Run from ${TB2026_ROOT}:"
    echo "  cmake -S k4SiWEcalReco -B k4SiWEcalReco/build && \\"
    echo "    cmake --build k4SiWEcalReco/build -j\$(nproc)"
    exit 1
fi

export LD_LIBRARY_PATH="${K4RECO_BUILD}:${LD_LIBRARY_PATH:-}"
export PYTHONPATH="${K4RECO_BUILD}/genConfDir:${TB2026_ROOT}:${REPO_ROOT}:${PYTHONPATH:-}"

# --------------------------------------------------------------------------- #
# Step 1: Convert simulation EDM4hep → ecal TTree
# --------------------------------------------------------------------------- #
echo ""
echo "=== Step 1/2: sim → ecal tree ==="
python3 -m analysis.sim_to_ecal_tree \
    --input  "${INPUT_EDM4HEP}" \
    --output "${ECAL_TREE}" \
    --collection "${COLLECTION}" \
    ${MAX_EVENTS} \
    --verbose

if [[ ! -f "${ECAL_TREE}" ]]; then
    echo "ERROR: converter did not produce ${ECAL_TREE}"
    exit 1
fi

# --------------------------------------------------------------------------- #
# Step 2: k4SiWEcalReco → EDM4hep / valtree with shower variables
# --------------------------------------------------------------------------- #
echo ""
echo "=== Step 2/2: k4SiWEcalReco (shower variables + event selection) ==="
python3 "${TB2026_ROOT}/k4SiWEcalReco/run_pid_batch.py" \
    --file   "${ECAL_TREE}" \
    --outdir "${OUT_DIR}" \
    --format "${FORMAT}" \
    ${VALIDATION_FLAG}

echo ""
echo "=== Done ==="
echo "Output directory: ${OUT_DIR}"
ls -lh "${OUT_DIR}"/ecal_sim*.root 2>/dev/null || true

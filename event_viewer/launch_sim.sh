#!/bin/bash
# ---------------------------------------------------------------------------
# launch_sim.sh
#
# Launch the event_viewer against simulation data produced by run_pid_sim.sh.
#
# Prerequisites
# -------------
#   a) run_pid_sim.sh completed successfully (ecal_sim.edm4hep.root exists)
#   b) key4hep environment sourced  (source init_key4hep.sh)
#   c) event_viewer Python dependencies installed:
#        bash setup_venv.sh
#
# Usage
# -----
#   cd /path/to/siwecal_k4sim
#   source init_key4hep.sh
#   bash event_viewer/launch_sim.sh [--port 8050] [--file PATH]
#
# Then open in a browser (or via SSH tunnel):
#   http://localhost:8050
# ---------------------------------------------------------------------------

set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

# Default: scan all pipeline outputs via sim_settings.yml (no --data-dir needed)
PORT=8050
FILE_ARG=""
DATA_DIR_ARG=""
DEBUG_FLAG=""

while [[ $# -gt 0 ]]; do
    case "$1" in
        --port)        PORT="$2"; shift 2 ;;
        --port=*)      PORT="${1#*=}"; shift ;;
        --file)        FILE_ARG="--file $2"; shift 2 ;;
        --file=*)      FILE_ARG="--file ${1#*=}"; shift ;;
        --data-dir)    DATA_DIR_ARG="--data-dir $2"; shift 2 ;;
        --debug)       DEBUG_FLAG="--debug"; shift ;;
        -h|--help)
            sed -n '/^# Usage/,/^# Then/p' "$0"
            exit 0 ;;
        *)
            echo "Unknown flag: $1"; exit 1 ;;
    esac
done

# If no --file given, auto-select the best available output from any pipeline
if [[ -z "${FILE_ARG}" ]]; then
    for candidate in \
        "${REPO_ROOT}/gaudi_jobs/1_e54_beam_pipeline/ecal_sim.edm4hep.root" \
        "${REPO_ROOT}/gaudi_jobs/1_mu_beam_pipeline/ecal_sim.edm4hep.root" \
        "${REPO_ROOT}/gaudi_jobs/1_e54_beam_pipeline/ecal_sim.valtree.root" \
        "${REPO_ROOT}/gaudi_jobs/1_mu_beam_pipeline/ecal_sim.valtree.root"; do
        if [[ -f "${candidate}" ]]; then
            FILE_ARG="--file ${candidate}"
            echo "[launch_sim] Auto-selected: ${candidate}"
            break
        fi
    done
    if [[ -z "${FILE_ARG}" ]]; then
        echo "WARNING: no k4SiWEcalReco output found in ${REPO_ROOT}/gaudi_jobs."
        echo "  Run  bash analysis/run_pid_sim.sh  first."
        echo "  Launching viewer without a pre-loaded file."
    fi
fi

# Activate the local virtualenv if it exists
VENV="${REPO_ROOT}/.venv-viewer"
if [[ -f "${VENV}/bin/activate" ]]; then
    # shellcheck disable=SC1090
    source "${VENV}/bin/activate"
fi

# Point the viewer at this repo's sim_settings.yml so data/geometry paths
# resolve against this repo root.
export SIWECAL_SETTINGS="${REPO_ROOT}/event_viewer/sim_settings.yml"

# Make the event_viewer package importable.
export PYTHONPATH="${REPO_ROOT}:${PYTHONPATH:-}"

echo "[launch_sim] Starting event_viewer on http://localhost:${PORT}"
echo "[launch_sim] Settings: ${SIWECAL_SETTINGS}"
echo ""

# Run from the siwecal_k4sim root so relative paths in scripts work naturally.
cd "${REPO_ROOT}"
python3 -m event_viewer \
    --port "${PORT}" \
    ${DATA_DIR_ARG} \
    ${FILE_ARG} \
    ${DEBUG_FLAG}

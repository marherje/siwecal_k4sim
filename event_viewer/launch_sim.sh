#!/bin/bash
# ---------------------------------------------------------------------------
# launch_sim.sh
#
# Launch the siwecal-tb2026 event_viewer against simulation data produced by
# run_pid_sim.sh.
#
# Prerequisites
# -------------
#   a) run_pid_sim.sh completed successfully (ecal_sim.edm4hep.root exists)
#   b) key4hep environment sourced  (source init_key4hep.sh)
#   c) event_viewer Python dependencies installed in siwecal-tb2026:
#        cd ../siwecal-tb2026
#        python -m venv .venv && source .venv/bin/activate
#        pip install -r requirements.txt
#
# Usage
# -----
#   cd /path/to/siwecal_k4sim
#   source init_key4hep.sh
#   bash event_viewer/launch_sim.sh [--port 8050] [--file ecal_sim.edm4hep.root]
#
# Then open in a browser (or via SSH tunnel):
#   http://localhost:8050
# ---------------------------------------------------------------------------

set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
TB2026_ROOT="$(cd "${REPO_ROOT}/../siwecal-tb2026" 2>/dev/null && pwd)" || {
    echo "ERROR: siwecal-tb2026 not found at ${REPO_ROOT}/../siwecal-tb2026"
    exit 1
}

DATA_DIR="${REPO_ROOT}/gaudi_jobs/1_mu_beam_pipeline"
PORT=8050
FILE_ARG=""
DEBUG_FLAG=""

while [[ $# -gt 0 ]]; do
    case "$1" in
        --port)        PORT="$2"; shift 2 ;;
        --port=*)      PORT="${1#*=}"; shift ;;
        --file)        FILE_ARG="--file $2"; shift 2 ;;
        --file=*)      FILE_ARG="--file ${1#*=}"; shift ;;
        --data-dir)    DATA_DIR="$2"; shift 2 ;;
        --debug)       DEBUG_FLAG="--debug"; shift ;;
        -h|--help)
            sed -n '/^# Usage/,/^# Then/p' "$0"
            exit 0 ;;
        *)
            echo "Unknown flag: $1"; exit 1 ;;
    esac
done

# If no --file given, auto-select the best available output
if [[ -z "${FILE_ARG}" ]]; then
    for candidate in \
        "${DATA_DIR}/ecal_sim.edm4hep.root" \
        "${DATA_DIR}/ecal_sim.valtree.root"; do
        if [[ -f "${candidate}" ]]; then
            FILE_ARG="--file ${candidate}"
            echo "[launch_sim] Auto-selected: ${candidate}"
            break
        fi
    done
    if [[ -z "${FILE_ARG}" ]]; then
        echo "WARNING: no k4SiWEcalReco output found in ${DATA_DIR}."
        echo "  Run  bash analysis/run_pid_sim.sh  first."
        echo "  Launching viewer without a pre-loaded file."
    fi
fi

# Activate the siwecal-tb2026 virtualenv if it exists
VENV="${TB2026_ROOT}/.venv"
if [[ -f "${VENV}/bin/activate" ]]; then
    # shellcheck disable=SC1090
    source "${VENV}/bin/activate"
fi

# Expose the geometry from this repo to the viewer
export SIWECAL_GEOMETRY_DIR="${REPO_ROOT}/geometry"

# Add both repos to PYTHONPATH so imports resolve regardless of cwd
export PYTHONPATH="${TB2026_ROOT}:${REPO_ROOT}:${PYTHONPATH:-}"

# Point siwecal_common to the simulation settings if the file exists
SIM_SETTINGS="${REPO_ROOT}/event_viewer/sim_settings.yml"
if [[ -f "${SIM_SETTINGS}" ]]; then
    export SIWECAL_SETTINGS="${SIM_SETTINGS}"
fi

echo "[launch_sim] Starting event_viewer on http://localhost:${PORT}"
echo "[launch_sim] Data dir: ${DATA_DIR}"
echo ""

cd "${TB2026_ROOT}"
python3 -m event_viewer \
    --data-dir "${DATA_DIR}" \
    --port "${PORT}" \
    ${FILE_ARG} \
    ${DEBUG_FLAG}

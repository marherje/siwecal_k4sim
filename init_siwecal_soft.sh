#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# init_siwecal_soft.sh
#
# Full environment setup for siwecal_k4sim: simulation, reconstruction,
# analysis, event_display (TEve) and event_viewer (Dash web app).
#
# Source this file — do NOT execute it:
#   source init_siwecal_soft.sh
#
# First-time setup (run once, not needed afterwards):
#   bash setup_venv.sh
# ---------------------------------------------------------------------------

REPO_ROOT="$( cd "$( dirname "${BASH_SOURCE[0]:-$0}" )" && pwd )"
TB2026_ROOT="$( cd "${REPO_ROOT}/../siwecal-tb2026" 2>/dev/null && pwd )" || TB2026_ROOT=""

# ---------------------------------------------------------------------------
# 1. key4hep stack (Gaudi, k4FWCore, DD4hep, ROOT, numpy, uproot, ...)
#    Release pinned in .key4hep-release (matching the siwecal-tb2026 checkout,
#    whose plugin is loaded in section 4); override per-shell with
#    KEY4HEP_RELEASE. Resolution is delegated to init_key4hep.sh so there is
#    only one copy of it in the repo.
# ---------------------------------------------------------------------------
source "${REPO_ROOT}/init_key4hep.sh" || return 1

# ---------------------------------------------------------------------------
# 2. Local build: SND_reco Gaudi algorithms + DD4hep detector plugin
# ---------------------------------------------------------------------------
export LD_LIBRARY_PATH="${REPO_ROOT}/install/lib64:${REPO_ROOT}/install/lib:${LD_LIBRARY_PATH}"
export PYTHONPATH="${REPO_ROOT}/install/lib64:${REPO_ROOT}/install/lib:${REPO_ROOT}/install/python:${REPO_ROOT}:${PYTHONPATH}"

# ---------------------------------------------------------------------------
# 3. Event viewer virtualenv (adds dash/plotly on top of key4hep)
# ---------------------------------------------------------------------------
VENV="${REPO_ROOT}/.venv-viewer"
if [[ -f "${VENV}/bin/activate" ]]; then
    source "${VENV}/bin/activate"
else
    echo "[siwecal] WARNING: event_viewer venv missing — run bash setup_venv.sh first (see README)"
fi

# ---------------------------------------------------------------------------
# 4. siwecal-tb2026: k4SiWEcalReco plugin + siwecal_validation
# ---------------------------------------------------------------------------
if [[ -n "${TB2026_ROOT}" ]]; then
    export PYTHONPATH="${TB2026_ROOT}:${PYTHONPATH}"

    # k4SiWEcalReco Gaudi plugin (shower variables: barz, moliere, ...)
    K4RECO_BUILD="${TB2026_ROOT}/gaudi_source/build"
    if [[ -f "${K4RECO_BUILD}/libk4SiWEcalRecoPlugins.so" ]]; then
        export LD_LIBRARY_PATH="${K4RECO_BUILD}:${LD_LIBRARY_PATH}"
        export PYTHONPATH="${K4RECO_BUILD}/genConfDir:${PYTHONPATH}"
    else
        echo "[siwecal] WARNING: k4SiWEcalReco not built — run bash setup_venv.sh first (see README)"
    fi
else
    echo "[siwecal] WARNING: siwecal-tb2026 not found at ${REPO_ROOT}/../siwecal-tb2026 (needed for k4SiWEcalReco)"
fi

# ---------------------------------------------------------------------------
# Summary
# ---------------------------------------------------------------------------
echo "[siwecal] Environment ready."
echo "[siwecal]   REPO_ROOT  = ${REPO_ROOT}"
[[ -n "${TB2026_ROOT}" ]] && echo "[siwecal]   TB2026_ROOT = ${TB2026_ROOT}"
echo "[siwecal]   k4run, ddsim, python3 -m event_viewer, python event_display/event_display_eve.py"

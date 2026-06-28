#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# setup_venv.sh  —  ONE-TIME setup (run once per machine/account)
#
# Creates the event_viewer virtualenv at .venv-viewer inside this repo.
#
# Usage:
#   source init_siwecal_soft.sh    # load key4hep first
#   bash setup_venv.sh
# ---------------------------------------------------------------------------
set -euo pipefail

REPO_ROOT="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

# Guard: key4hep must be active so the venv is created with Python 3.13,
# not the system Python 3.9 which lacks podio/uproot.
PY_MINOR=$(python3 -c "import sys; print(sys.version_info.minor)")
if [[ "$(python3 -c 'import sys; print(sys.version_info.major)')" -lt 3 || "${PY_MINOR}" -lt 13 ]]; then
    echo "ERROR: Python 3.13+ required (got $(python3 --version))."
    echo "       Source key4hep first:  source init_siwecal_soft.sh"
    exit 1
fi

echo "=== Creating event_viewer virtualenv ==="
VENV="${REPO_ROOT}/.venv-viewer"
if [[ -f "${VENV}/bin/activate" ]]; then
    # Check the venv also uses Python 3.13 (in case it was created with 3.9 before)
    VENV_PY_MINOR=$("${VENV}/bin/python3" -c "import sys; print(sys.version_info.minor)" 2>/dev/null || echo "0")
    if [[ "${VENV_PY_MINOR}" -lt 13 ]]; then
        echo "  Existing venv uses Python < 3.13 — recreating with Python 3.13..."
        rm -rf "${VENV}"
    else
        echo "  Venv already exists (Python 3.${VENV_PY_MINOR}) — skipping creation."
    fi
fi
if [[ ! -f "${VENV}/bin/activate" ]]; then
    python3 -m venv --system-site-packages "${VENV}"
    echo "  Created: ${VENV} ($(python3 --version))"
fi

source "${VENV}/bin/activate"
pip install --quiet -r "${REPO_ROOT}/event_viewer/requirements.txt"
python3 -c "import dash, plotly; print(f'  dash {dash.__version__}, plotly {plotly.__version__} — OK')"
deactivate

echo ""
echo "=== Setup complete ==="
echo "From now on, just run:  source init_siwecal_soft.sh"

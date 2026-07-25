# Source the key4hep stack. Single source of truth for the release is
# .key4hep-release (shared with build.sh / init_siwecal_soft.sh and matching the
# siwecal-tb2026 checkout, whose k4SiWEcalReco plugin this project loads);
# override per-shell with the KEY4HEP_RELEASE env var.
_K4_ROOT="$( cd "$( dirname "${BASH_SOURCE[0]:-$0}" )" && pwd )"
KEY4HEP_RELEASE="${KEY4HEP_RELEASE:-$(cat "${_K4_ROOT}/.key4hep-release" 2>/dev/null || echo 2026-04-08)}"
source /cvmfs/sw.hsf.org/key4hep/setup.sh -r "${KEY4HEP_RELEASE}"

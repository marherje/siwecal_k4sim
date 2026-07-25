# Source the key4hep stack. Single source of truth for the release is
# .key4hep-release (versioned in the repo, matching the siwecal-tb2026 checkout
# whose k4SiWEcalReco plugin this project loads); override per-shell with the
# KEY4HEP_RELEASE env var.
#
# build.sh and init_siwecal_soft.sh source THIS file instead of repeating the
# resolution: three copies of it existed, each with its own hardcoded fallback,
# so updating the pin in one place left the others silently on the old release.
_K4_ROOT="$( cd "$( dirname "${BASH_SOURCE[0]:-$0}" )" && pwd )"
_K4_PIN="${_K4_ROOT}/.key4hep-release"
if [[ -z "${KEY4HEP_RELEASE:-}" ]]; then
    if [[ ! -r "${_K4_PIN}" ]]; then
        # Deliberately no fallback release here: the previous "|| echo <date>"
        # meant a missing or unreadable pin started an unintended stack without
        # a word, which is indistinguishable from a healthy setup until jobs
        # disagree with each other.
        echo "ERROR: missing key4hep pin ${_K4_PIN} — it is part of the repo." >&2
        echo "       Restore it (git checkout -- .key4hep-release) or set KEY4HEP_RELEASE." >&2
        return 1 2>/dev/null || exit 1
    fi
    KEY4HEP_RELEASE="$(tr -d '[:space:]' < "${_K4_PIN}")"
fi
source /cvmfs/sw.hsf.org/key4hep/setup.sh -r "${KEY4HEP_RELEASE}"

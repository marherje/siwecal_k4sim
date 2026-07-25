#!/bin/bash
# Shared processing pipeline for the 2026 PID campaign (e-/mu-/pi-, 74 & 99 GeV,
# beam centre (1.139, 1.164) mm, 50k events per point).
#
# The condor chunk jobs (simulation/run_script/launch_beam_pid2026.sh) already
# did ddsim + digitisation + ecal-tree conversion and left one small ecal TTree
# per chunk in  Processed/chunks/.  This script therefore only has to:
#
#   1. hadd the chunk trees of the requested sample  ->  ecal_<label>.root
#      (for --energy merged: the chunks of BOTH 74 and 99 GeV, which is exactly
#       the "merge at the end" strategy — no reprocessing, no podio metadata
#       games; the per-chunk 'run' branch keeps energy+chunk identity)
#   2. run k4SiWEcalReco (run_pid_batch.py) -> shower variables
#   3. stage the products to  Processed/  as <label>_ecal.{root,edm4hep.root,valtree.root}
#
# Work happens under WORKDIR (default /tmp/$USER/...), never on AFS: the merged
# trees are ~0.8-1.6 GB and the AFS home volume has ~1.5 GB free.
#
# Usage:
#   bash pipeline_common.sh --particle e- --energy 74
#   bash pipeline_common.sh --particle pi- --energy merged [--allow-partial]

set -uo pipefail

PARTICLE=""
ENERGY=""
ALLOW_PARTIAL=0
NCHUNKS=${NCHUNKS:-25}
FORMAT=${FORMAT:-both}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --particle)      PARTICLE="$2"; shift 2 ;;
        --energy)        ENERGY="$2"; shift 2 ;;
        --allow-partial) ALLOW_PARTIAL=1; shift ;;
        --nchunks)       NCHUNKS="$2"; shift 2 ;;
        --format)        FORMAT="$2"; shift 2 ;;
        -h|--help)       sed -n '2,26p' "$0"; exit 0 ;;
        *) echo "Unknown flag: $1"; exit 1 ;;
    esac
done

[[ -n "${PARTICLE}" && -n "${ENERGY}" ]] || { echo "ERROR: --particle and --energy are required"; exit 1; }

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
TB2026_ROOT="$(cd "${REPO_ROOT}/../siwecal-tb2026" 2>/dev/null && pwd)" || {
    echo "ERROR: siwecal-tb2026 not found next to ${REPO_ROOT}"; exit 1; }

EOS_BASE="${EOS_BASE:-/eos/experiment/drdcalo/siw-ecal/TB2026-06/Simulation}"
CHUNK_DIR="${EOS_BASE}/Processed/chunks"
PROCESSED="${EOS_BASE}/Processed"
WORKDIR="${WORKDIR:-/tmp/${USER}/siwecal_pid2026}"
mkdir -p "${WORKDIR}"

POS_X=1.139
POS_Y=1.164
SIGMA_E=0.02
case "${PARTICLE}" in
    "e-")  SIGX=13.75; SIGY=8.25  ;;
    "mu-") SIGX=38.5;  SIGY=46.75 ;;
    "pi-") SIGX=24.75; SIGY=13.75 ;;
    *) echo "ERROR: unknown particle '${PARTICLE}'"; exit 1 ;;
esac

point_label() {  # $1 = energy in GeV
    echo "beam_${PARTICLE}_${1}GeV_xy_${POS_X}_${POS_Y}_sigx${SIGX}_sigy${SIGY}_sigE${SIGMA_E}"
}

if [[ "${ENERGY}" == "merged" ]]; then
    ENERGIES="74 99"
    LABEL="beam_${PARTICLE}_74-99GeV_xy_${POS_X}_${POS_Y}_sigx${SIGX}_sigy${SIGY}_sigE${SIGMA_E}"
else
    ENERGIES="${ENERGY}"
    LABEL="$(point_label "${ENERGY}")"
fi

echo "=== PID 2026 pipeline ==="
echo "  particle : ${PARTICLE}"
echo "  energy   : ${ENERGY} GeV"
echo "  label    : ${LABEL}"
echo "  chunks   : ${CHUNK_DIR}"
echo "  workdir  : ${WORKDIR}"

# --------------------------------------------------------------------------- #
# Step 1: collect + hadd the chunk trees
# --------------------------------------------------------------------------- #
CHUNK_FILES=()
for e in ${ENERGIES}; do
    pl="$(point_label "${e}")"
    found=()
    while IFS= read -r f; do [[ -n "$f" ]] && found+=("$f"); done \
        < <(ls -1 "${CHUNK_DIR}/${pl}_c"*"_ecal.root" 2>/dev/null | sort)
    echo "  ${e} GeV: ${#found[@]} / ${NCHUNKS} chunk trees"
    if (( ${#found[@]} != NCHUNKS )) && (( ALLOW_PARTIAL == 0 )); then
        echo "ERROR: expected ${NCHUNKS} chunk trees for ${pl}, found ${#found[@]}."
        echo "       Wait for the condor jobs, resubmit the missing chunks, or pass --allow-partial."
        exit 1
    fi
    (( ${#found[@]} == 0 )) && { echo "ERROR: no chunk trees at all for ${pl}"; exit 1; }
    CHUNK_FILES+=("${found[@]}")
done

MERGED_TREE="${WORKDIR}/ecal_${LABEL}.root"
echo ""
echo "=== Step 1/3: hadd ${#CHUNK_FILES[@]} chunk trees -> $(basename "${MERGED_TREE}") ==="
hadd -f -k "${MERGED_TREE}" "${CHUNK_FILES[@]}" > "${WORKDIR}/hadd_${LABEL}.log" 2>&1 || {
    echo "ERROR: hadd failed — see ${WORKDIR}/hadd_${LABEL}.log"; exit 1; }
python3 -c "
import ROOT, sys
f = ROOT.TFile.Open('${MERGED_TREE}')
t = f.Get('ecal')
print(f'  merged entries: {t.GetEntries()}')
sys.exit(0 if t.GetEntries() > 0 else 1)
" || { echo "ERROR: merged tree is empty"; exit 1; }

# --------------------------------------------------------------------------- #
# Step 2: k4SiWEcalReco (shower variables + event selection)
# --------------------------------------------------------------------------- #
K4RECO_BUILD="${TB2026_ROOT}/gaudi_source/build"
if [[ ! -f "${K4RECO_BUILD}/libk4SiWEcalRecoPlugins.so" ]]; then
    echo "ERROR: k4SiWEcalReco not built. From ${TB2026_ROOT}:"
    echo "  cmake -S gaudi_source -B gaudi_source/build && cmake --build gaudi_source/build -j\$(nproc)"
    exit 1
fi
export LD_LIBRARY_PATH="${K4RECO_BUILD}:${LD_LIBRARY_PATH:-}"
export PYTHONPATH="${K4RECO_BUILD}/genConfDir:${TB2026_ROOT}:${REPO_ROOT}:${PYTHONPATH:-}"

echo ""
echo "=== Step 2/3: k4SiWEcalReco (format=${FORMAT}) ==="
python3 "${TB2026_ROOT}/gaudi_jobs/run_pid_batch.py" \
    --file   "${MERGED_TREE}" \
    --outdir "${WORKDIR}" \
    --format "${FORMAT}" 2>&1 | grep -v "^TCling::LoadPCM"
RECO_RC=${PIPESTATUS[0]}
if [[ ${RECO_RC} -ne 0 ]]; then
    echo "ERROR: run_pid_batch.py failed (exit ${RECO_RC})"; exit 1
fi

# --------------------------------------------------------------------------- #
# Step 3: stage products to Processed/ (naming as in the 1_e74 pipeline)
# --------------------------------------------------------------------------- #
echo ""
echo "=== Step 3/3: staging to ${PROCESSED}/ ==="
staged=0
for spec in "root:${MERGED_TREE}" \
            "edm4hep.root:${WORKDIR}/ecal_${LABEL}.edm4hep.root" \
            "valtree.root:${WORKDIR}/ecal_${LABEL}.valtree.root"; do
    ext="${spec%%:*}"; src="${spec#*:}"
    [[ -f "${src}" ]] || continue
    dst="${PROCESSED}/${LABEL}_ecal.${ext}"
    cp -f "${src}" "${dst}" && rm -f "${src}" && staged=$(( staged + 1 )) \
        || { echo "ERROR: could not stage ${src} -> ${dst}"; exit 1; }
    echo "  ${dst}"
done
rm -f "${WORKDIR}/.${LABEL}.pid.tmp.root"   # run_pid_batch strips the 'ecal_' prefix

echo ""
echo "=== Done — ${staged} file(s) in ${PROCESSED}/ ==="
ls -lh "${PROCESSED}/${LABEL}_ecal."* 2>/dev/null || true

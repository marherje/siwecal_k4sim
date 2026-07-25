#!/bin/bash
# Run the three PID-2026 pipelines of one particle (74, 99, merged 74+99) and,
# only if all three produced their final samples, retire that particle's chunk
# trees from Processed/chunks/.
#
# Why the cleanup lives here and not in pipeline_common.sh: the 74 GeV chunks
# are consumed twice — by the 74 GeV pipeline and by the merged one — so no
# single pipeline may delete them. This driver is the only place that knows all
# consumers are done.
#
# Only the ecal chunk TREES are removed. The raw simulation chunks in
# Generated/chunks/ are kept: they are the reprocessable source, and
# regenerating them costs ~26 h of CPU per point.
#
# Usage:
#   bash run_particle.sh --particle e- [--keep-chunks] [--allow-partial]

set -uo pipefail

PARTICLE=""
KEEP_CHUNKS=0
PASSTHRU=()

while [[ $# -gt 0 ]]; do
    case "$1" in
        --particle)    PARTICLE="$2"; shift 2 ;;
        --keep-chunks) KEEP_CHUNKS=1; shift ;;
        *)             PASSTHRU+=("$1"); shift ;;
    esac
done
[[ -n "${PARTICLE}" ]] || { echo "ERROR: --particle is required"; exit 1; }

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
EOS_BASE="${EOS_BASE:-/eos/experiment/drdcalo/siw-ecal/TB2026-06/Simulation}"
CHUNK_DIR="${EOS_BASE}/Processed/chunks"
PROCESSED="${EOS_BASE}/Processed"

POS_X=1.139; POS_Y=1.164; SIGMA_E=0.02
case "${PARTICLE}" in
    "e-")  SIGX=13.75; SIGY=8.25  ;;
    "mu-") SIGX=38.5;  SIGY=46.75 ;;
    "pi-") SIGX=24.75; SIGY=13.75 ;;
    *) echo "ERROR: unknown particle '${PARTICLE}'"; exit 1 ;;
esac
SUFFIX="xy_${POS_X}_${POS_Y}_sigx${SIGX}_sigy${SIGY}_sigE${SIGMA_E}"

failed=0
for tag in 74 99 merged; do
    echo ""
    echo "################ ${PARTICLE} — ${tag} ################"
    bash "${HERE}/pipeline_common.sh" --particle "${PARTICLE}" --energy "${tag}" \
         ${PASSTHRU[@]+"${PASSTHRU[@]}"}
    rc=$?
    if [[ ${rc} -ne 0 ]]; then
        echo "!! pipeline ${PARTICLE} ${tag} failed (exit ${rc})"
        failed=$(( failed + 1 ))
    fi
done

if (( failed > 0 )); then
    echo ""
    echo "=== ${failed}/3 pipeline(s) failed for ${PARTICLE} — chunk trees kept for retry ==="
    exit 1
fi

# --------------------------------------------------------------------------- #
# Retire the chunk trees, but only against evidence that all three samples
# really landed in Processed/ with a non-trivial size.
# --------------------------------------------------------------------------- #
for tag in 74 99 74-99; do
    out="${PROCESSED}/beam_${PARTICLE}_${tag}GeV_${SUFFIX}_ecal.valtree.root"
    if [[ ! -s "${out}" ]]; then
        echo "=== Final sample missing (${out}) — chunk trees kept ==="
        exit 1
    fi
done

if (( KEEP_CHUNKS == 1 )); then
    echo ""
    echo "=== --keep-chunks: leaving ${PARTICLE} chunk trees in ${CHUNK_DIR} ==="
    exit 0
fi

echo ""
echo "=== All three ${PARTICLE} samples verified — retiring chunk trees ==="
removed=0
for e in 74 99; do
    for f in "${CHUNK_DIR}/beam_${PARTICLE}_${e}GeV_${SUFFIX}"_c*_ecal.root; do
        [[ -e "${f}" ]] || continue
        rm -f "${f}" && removed=$(( removed + 1 ))
    done
done
echo "  removed ${removed} chunk tree(s) from ${CHUNK_DIR}"
echo "  raw sim chunks kept in ${EOS_BASE}/Generated/chunks/ for reprocessing"
echo ""
ls -lh "${PROCESSED}/beam_${PARTICLE}_"*"GeV_${SUFFIX}_ecal."* 2>/dev/null

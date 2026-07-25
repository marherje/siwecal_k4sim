#!/bin/bash
# Delete the logs of PID-2026 chunk jobs that are finished, to keep the AFS home
# volume from filling while a campaign runs (each electron chunk leaves ~700 KB
# of ddsim log, and 150 chunks + reprocess jobs add up).
#
# A chunk counts as finished when its raw sim is on EOS *and* no job for it is
# left in the condor queue. That deliberately keeps:
#   * logs of running/idle/held jobs — still being written, or needed to
#     diagnose the hold
#   * logs of chunks that left no simulation behind — the only record of why
#
# Steer files (.mac/.py/.sh/.sub) are never touched: the condor shell script
# reads the GPS macro and steering file from steer/ at run time, so removing
# them would break jobs that are still queued.
#
# Usage:
#   ./prune_logs.sh            # prune
#   ./prune_logs.sh --dry-run  # list what would go, delete nothing

set -uo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")"

DRY=0
[[ "${1:-}" == "--dry-run" ]] && DRY=1

SIM_DIR=/eos/experiment/drdcalo/siw-ecal/TB2026-06/Simulation/Generated/chunks
LOG=log

# Labels that still have a job in the queue (any status, including held).
queued=$(condor_q "$USER" -af Cmd 2>/dev/null \
         | grep -oE '(runddsim_QGSP_BERT_|reproc_)[^/]+\.sh' \
         | sed -E 's/^(runddsim_QGSP_BERT_|reproc_)//; s/\.sh$//' | sort -u)

freed=0
pruned=0
kept_active=0
kept_nosim=0

for f in ${LOG}/beam_*_c[0-9][0-9][0-9].log; do
    [[ -e "$f" ]] || continue
    label=$(basename "$f" .log)

    if grep -qxF "${label}" <<< "${queued}"; then
        kept_active=$(( kept_active + 1 )); continue
    fi
    if [[ ! -s "${SIM_DIR}/output_${label}.edm4hep.root" ]]; then
        kept_nosim=$(( kept_nosim + 1 )); continue
    fi

    for victim in "${LOG}/${label}.log" \
                  "${LOG}/reproc_${label}.log" \
                  "${LOG}/runddsim_QGSP_BERT_${label}.log" \
                  "${LOG}/reproc_${label}.condor.log" \
                  "${LOG}/outfile_runddsim_QGSP_BERT_${label}.txt" \
                  "${LOG}/errors_runddsim_QGSP_BERT_${label}.txt" \
                  "${LOG}/outfile_reproc_${label}.txt" \
                  "${LOG}/errors_reproc_${label}.txt"; do
        [[ -e "${victim}" ]] || continue
        sz=$(stat -c%s "${victim}" 2>/dev/null || echo 0)
        if (( DRY == 1 )); then
            echo "  would remove $(basename "${victim}") ($(( sz / 1024 )) KB)"
        else
            rm -f "${victim}" && freed=$(( freed + sz ))
        fi
    done
    pruned=$(( pruned + 1 ))
done

if (( DRY == 1 )); then
    echo "DRY RUN: ${pruned} finished chunk(s) would be pruned; ${kept_active} active, ${kept_nosim} without sim kept"
else
    echo "PRUNED ${pruned} finished chunk(s), freed $(( freed / 1024 / 1024 )) MB; kept ${kept_active} active + ${kept_nosim} without sim"
fi

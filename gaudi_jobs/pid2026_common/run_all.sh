#!/bin/bash
# Process the whole PID-2026 campaign: for each particle run the 74, 99 and
# merged pipelines and retire that particle's chunk trees when all three are
# verified in Processed/ (see run_particle.sh).
#
# Usage:
#   source init_siwecal_soft.sh
#   bash gaudi_jobs/pid2026_common/run_all.sh [--keep-chunks] [--allow-partial]
#   PARTICLES="e- pi-" bash gaudi_jobs/pid2026_common/run_all.sh

set -uo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PARTICLES=${PARTICLES:-"e- mu- pi-"}

rc_total=0
for p in ${PARTICLES}; do
    bash "${HERE}/run_particle.sh" --particle "${p}" "$@" || rc_total=1
done

echo ""
if (( rc_total == 0 )); then
    echo "=== Campaign processed: all requested particles done ==="
else
    echo "=== Some particles failed — their chunk trees were kept for retry ==="
fi
exit ${rc_total}

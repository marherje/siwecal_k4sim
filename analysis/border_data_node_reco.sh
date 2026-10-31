#!/bin/bash
# Same as border_data_node.sh, but with the Reconstruction campaign's own recipe (generate_reco_dag.py reco.sh):
# k4run options/run_pid.py (ECAL_HIT_MIP_CUT 0.5, default mip05_/mip1_ blocks), then pid_to_valtree.py. Unlike
# run_pid_batch.py, it keeps the events without hits, so the valtree stays row-aligned with ecal_<run>.root. Used for
# the runs that have empty events (107, 111, 116, 134).
#
#   bash analysis/border_data_node_reco.sh <run> <reco_base> <out_base>
RUN=${1:?run}; RECO=${2:?reco_base}; O=${3:?out_base}
R="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
T=${TB_REPO:-$(cd "$R/../siwecal-tb2026" && pwd)}
W=${WORK:-${TMPDIR:-/tmp}/border_data_reco_$USER}/$RUN
IN=$RECO/$RUN/ecal_$RUN.root

[ -f "$IN" ] || { echo "missing input $IN" >&2; exit 1; }
mkdir -p "$W"
set -o pipefail
(cd "$T" && source "$T/setup.sh" >/dev/null 2>&1 && \
  ECAL_FILE="$IN" ECAL_PID_OUT="$W/ecal_$RUN.edm4hep.root" ECAL_HIT_MIP_CUT=0.5 \
  k4run gaudi_source/options/run_pid.py > "$W/pid.log" 2>&1 && \
  python3 gaudi_jobs/pid_to_valtree.py --force "$W/ecal_$RUN.edm4hep.root" >> "$W/pid.log" 2>&1) \
  || { echo "pid FAILED $RUN"; tail -40 "$W/pid.log"; exit 1; }
[ -d "$O/$RUN" ] || mkdir "$O/$RUN"
for kind in edm4hep valtree; do
  f=$W/ecal_$RUN.$kind.root
  [ -s "$f" ] || { echo "no $kind for $RUN"; tail -5 "$W/pid.log"; exit 1; }
  cp "$f" "$O/$RUN/" || { echo "copy FAILED $RUN $kind"; exit 1; }
  [ "$(stat -c %s "$f")" = "$(stat -c %s "$O/$RUN/$(basename "$f")")" ] || { echo "size mismatch $RUN $kind"; exit 1; }
done
cp "$W/pid.log" "$O/$RUN/"
rm -rf "$W"
echo "DONE $RUN"

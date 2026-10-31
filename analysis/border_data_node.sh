#!/bin/bash
# One data run (node of analysis/make_border_data_dag.py): re-run the PID stage on Reconstruction/<run>/ecal_<run>.root
# with the current siwecal-tb2026 build, so that the EDM4hep and the valtree carry barxy_min_border_distance. Same
# settings as the Reconstruction campaign (run_pid.py: hit MIP cut 0.5, mip05_/mip1_ variant blocks), so every older
# branch is reproduced. The input is only read; the outputs go to a staging directory, from which they replace
# Reconstruction/<run>/ecal_<run>.{edm4hep,valtree}.root.
#
#   bash analysis/border_data_node.sh <run> <reco_base> <out_base>
#   in : <reco_base>/<run>/ecal_<run>.root
#   out: <out_base>/<run>/ecal_<run>.{edm4hep,valtree}.root (+ pid.log)
RUN=${1:?run}; RECO=${2:?reco_base}; O=${3:?out_base}
R="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
T=${TB_REPO:-$(cd "$R/../siwecal-tb2026" && pwd)}
W=${WORK:-${TMPDIR:-/tmp}/border_data_$USER}/$RUN
IN=$RECO/$RUN/ecal_$RUN.root

[ -f "$IN" ] || { echo "missing input $IN" >&2; exit 1; }
mkdir -p "$W"
set -o pipefail
(cd "$T" && source "$T/setup.sh" >/dev/null 2>&1 && \
  python3 gaudi_jobs/run_pid_batch.py --file "$IN" --outdir "$W" --format both --validation --hit-mip-cut 0.5 \
  > "$W/pid.log" 2>&1) || { echo "pid FAILED $RUN"; tail -40 "$W/pid.log"; free -m; exit 1; }
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

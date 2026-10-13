#!/bin/bash
# One beam point after all its chunks (node of analysis/make_sim_dag.py): hadd the chunk ecal trees written by
# simulation/run_script/sim_digi_chunk.sh -> trees/ecal_<name>{,_simple}.root -> ECalPid edm4hep (hit MIP cut 0.5)
# -> valtree, for the digitised (digi) and undigitised (nodigi) samples; the chunk trees are removed at the end.
# Same PID/valtree steps and output layout as run_sim_chain_final.sh (final_v6).
#
#   bash analysis/sim_reco_group.sh <name> <out_base> <nchunks>
#   in : <out_base>/chunks/<name>/ecal_<name>_cNNN{,_simple}.root
#   out: <out_base>/trees/ecal_<name>{,_simple}.root, <out_base>/pid/<name>_{digi,nodigi}/*.{edm4hep,valtree}.root
NAME=${1:?name}; O=${2:?out_base}; NCH=${3:?nchunks}
R="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
T=${TB_REPO:-$(cd "$R/../siwecal-tb2026" && pwd)}
W=${WORK:-${TMPDIR:-/tmp}/sim_reco_$USER}/$NAME
C=$O/chunks/$NAME

source "$R/init_key4hep.sh" >/dev/null 2>&1
command -v hadd >/dev/null || { echo "hadd not on PATH" >&2; exit 1; }
set -o pipefail
for d in "$O/trees" "$O/pid" "$O/pid/${NAME}_digi" "$O/pid/${NAME}_nodigi"; do [ -d "$d" ] || mkdir "$d"; done
mkdir -p "$W"

for mode in real simple; do
  if [ $mode = simple ]; then sfx=_simple; out=$O/trees/ecal_${NAME}_simple.root; pid=$O/pid/${NAME}_nodigi
  else sfx=""; out=$O/trees/ecal_${NAME}.root; pid=$O/pid/${NAME}_digi; fi
  shopt -s nullglob; trees=("$C"/ecal_${NAME}_c[0-9][0-9][0-9]${sfx}.root); shopt -u nullglob
  [ ${#trees[@]} -eq "$NCH" ] || { echo "$NAME $mode: ${#trees[@]}/$NCH chunk trees in $C" >&2; exit 1; }
  hadd -f "$W/ecal.root" "${trees[@]}" > "$W/hadd_$mode.log" 2>&1 || { echo "hadd FAILED $NAME $mode"; tail -5 "$W/hadd_$mode.log"; exit 1; }
  cp "$W/ecal.root" "$out" && rm -f "$W/ecal.root" || { echo "copy FAILED $out"; exit 1; }
  (cd "$T" && source "$T/setup.sh" >/dev/null 2>&1 && \
    python3 gaudi_jobs/run_pid_batch.py --file "$out" --outdir "$W/pid_$mode" --format edm4hep --hit-mip-cut 0.5 \
    > "$W/pid_$mode.log" 2>&1) || { echo "pid FAILED $NAME $mode"; tail -5 "$W/pid_$mode.log"; exit 1; }
  rm -f "$pid"/*.root; cp "$W/pid_$mode"/*.edm4hep.root "$pid"/ || { echo "copy FAILED $pid"; exit 1; }
  (cd "$T" && source "$T/setup.sh" >/dev/null 2>&1 && \
    python3 gaudi_jobs/pid_to_valtree.py --force "$pid"/*.edm4hep.root > "$W/valtree_$mode.log" 2>&1) \
    || { echo "valtree FAILED $NAME $mode"; tail -5 "$W/valtree_$mode.log"; exit 1; }
done
rm -rf "$W" "$C"
echo "DONE $NAME"

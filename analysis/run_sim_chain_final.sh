#!/bin/bash
# One final_v6 group through the analysis chain: its ddsim chunks (simulation/run_script/launch_beam_final.sh)
# -> digitised (real, reconstructed with the th210 table) and undigitised (simple) ecal trees, one job3 per chunk,
# hadd'ed per mode -> ECalPid edm4hep (hit MIP cut 0.5) -> valtree. Same chain and defaults as run_sim_chain.sh.
#
#   bash analysis/run_sim_chain_final.sh <set> <pos> <E_GeV> [<out_base>]     e.g.  th230 E1 74
#   (condor: analysis/submit_sim_chain_final.sh; 25 chunks x 2 modes of job3 + PID on 50k events: up to ~8 h at 99 GeV)
#
# Input : Generated/final_v6/output_beam_e-_<E>GeV_<pos>_xy_*_c*.edm4hep.root (every chunk of the group)
# Output: <out_base>/trees/ecal_<set>_<pos>_e<E>{,_simple}.root
#         <out_base>/pid/<set>_<pos>_e<E>_{digi,nodigi}/*.{edm4hep,valtree}.root
# out_base defaults to Simulation/Processed/adc_vs_tb/final_v6. Group names match final_panels.py with
# FINAL_GROUPS=analysis/final_v6_groups.txt FINAL_SIM=<out_base>.
#
# Environment: DIGI_CALIB_DIR, RECO_TABLE, TB_REPO, PED_NOISE, WORK, FORCE, FORCE_DIGI as in run_sim_chain.sh;
#   IN_DIR (ddsim chunks, default Generated/final_v6), PARTICLE (e-), MIN_CHUNKS (refuse fewer, default 1),
#   NPAR (job3 runs in parallel, default 1).
SET=${1:?set}; POS=${2:?pos}; E=${3:?energy}
R="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
T=${TB_REPO:-$(cd "$R/../siwecal-tb2026" && pwd)}
B=/eos/experiment/drdcalo/siw-ecal/TB2026-06/Simulation
O=${4:-$B/Processed/adc_vs_tb/final_v6}
IN_DIR=${IN_DIR:-$B/Generated/final_v6}
export DIGI_CALIB_DIR=${DIGI_CALIB_DIR:-/eos/experiment/drdcalo/siw-ecal/TB2026-06/Calib/MuonCalib}
RECO_TABLE=${RECO_TABLE:-th210}
export PED_NOISE=${PED_NOISE:-1}
P=${PARTICLE:-e-}
case "$P" in e-) tag=${SET}_${POS}_e${E} ;; *) tag=${SET}_${POS}_${P%[-+]}${E} ;; esac
W=${WORK:-${TMPDIR:-/tmp}/sim_chain_$USER}/$tag

# key4hep before any `set -u` (it would leave k4run off the PATH), then the repo's own components.
source "$R/init_key4hep.sh" >/dev/null 2>&1
command -v k4run >/dev/null || { echo "k4run not on PATH" >&2; exit 1; }
set -o pipefail
export LD_LIBRARY_PATH=$R/install/lib64:$R/install/lib:$LD_LIBRARY_PATH
export PYTHONPATH=$R/install/lib64:$R/install/lib:$R/install/python:$PYTHONPATH

shopt -s nullglob
chunks=("$IN_DIR"/output_beam_${P}_${E}GeV_${POS}_xy_*_c[0-9][0-9][0-9].edm4hep.root)
shopt -u nullglob
[ ${#chunks[@]} -ge "${MIN_CHUNKS:-1}" ] || { echo "only ${#chunks[@]} ddsim chunks for $tag in $IN_DIR" >&2; exit 1; }
# one sample per group: every chunk must carry the same point name
pts=$(printf "%s\n" "${chunks[@]##*/}" | sed 's/_c[0-9]*\.edm4hep\.root$//' | sort -u)
[ "$(echo "$pts" | wc -l)" = 1 ] || { echo "chunks of different samples for $tag:"$'\n'"$pts" >&2; exit 1; }
echo "$tag: ${#chunks[@]} chunks of $pts"

[ -d "$DIGI_CALIB_DIR/mips/$RECO_TABLE" ] || { echo "no $RECO_TABLE table in $DIGI_CALIB_DIR" >&2; exit 1; }
# plain mkdir per level: mkdir -p walks up into EOS parents it cannot touch
for d in "$O" "$O/trees" "$O/pid" "$O/pid/${tag}_digi" "$O/pid/${tag}_nodigi"; do [ -d "$d" ] || mkdir "$d"; done
mkdir -p "$W"

for mode in real simple; do
  if [ $mode = simple ]; then
    out=$O/trees/ecal_${tag}_simple.root; pid=$O/pid/${tag}_nodigi; reco=""
  else
    out=$O/trees/ecal_${tag}.root; pid=$O/pid/${tag}_digi; reco=$RECO_TABLE
  fi
  redo=${FORCE:-0}; [ $mode = real ] && [ "${FORCE_DIGI:-0}" = 1 ] && redo=1
  if [ "$redo" = 1 ] || [ ! -s "$out" ]; then
    # one job3 per chunk, each in its own directory (job3 writes digitized.edm4hep.root in its cwd);
    # NPAR of them at a time (default 1, the condor jobs)
    digi_chunk() {  # chunk file -> $W/<mode>_<cNNN>/ecal.root
      local f=$1 c d; c=$(basename "$f" .edm4hep.root); c=${c##*_}; d="$W/${mode}_$c"; mkdir -p "$d"
      # the ecal tree is written by EcalTreeWriter in the same k4run (ECAL_TREE_OUTPUT)
      (cd "$d" && CALIB_THRESHOLD=$SET DIGI_MODE=$mode RECO_TABLE_THRESHOLD=$reco INPUT_FILE=$f \
        ECAL_TREE_OUTPUT=ecal.root k4run "$R/gaudi_jobs/pid2026_common/job3_digitize.py" > job3.log 2>&1) \
        || { echo "job3 FAILED $tag $mode $c"; tail -5 "$d/job3.log"; return 1; }
      rm -f "$d/digitized.edm4hep.root"
    }
    export -f digi_chunk; export W R SET mode reco tag
    printf '%s\n' "${chunks[@]}" | xargs -P "${NPAR:-1}" -I{} bash -c 'digi_chunk "$1"' _ {} || exit 1
    trees=()
    for f in "${chunks[@]}"; do c=$(basename "$f" .edm4hep.root); trees+=("$W/${mode}_${c##*_}/ecal.root"); done
    hadd -f "$W/ecal.root" "${trees[@]}" > "$W/hadd_$mode.log" 2>&1 || { echo "hadd FAILED $tag $mode"; exit 1; }
    cp "$W/ecal.root" "$out" && rm -f "$W/ecal.root" && rm -rf "$W/${mode}"_c[0-9]*
  fi
  if [ "$redo" = 1 ] || ! ls "$pid"/*.edm4hep.root >/dev/null 2>&1; then
    (cd "$T" && source "$T/setup.sh" >/dev/null 2>&1 && \
      python3 gaudi_jobs/run_pid_batch.py --file "$out" --outdir "$W/pid_$mode" --format edm4hep --hit-mip-cut 0.5 \
      > "$W/pid_$mode.log" 2>&1) || { echo "pid FAILED $tag $mode"; tail -5 "$W/pid_$mode.log"; exit 1; }
    cp "$W/pid_$mode"/*.edm4hep.root "$pid"/
  fi
  vforce=""; [ "$redo" = 1 ] && vforce=--force
  (cd "$T" && source "$T/setup.sh" >/dev/null 2>&1 && \
    python3 gaudi_jobs/pid_to_valtree.py $vforce "$pid"/*.edm4hep.root > "$W/valtree_$mode.log" 2>&1) \
    || { echo "valtree FAILED $tag $mode"; tail -5 "$W/valtree_$mode.log"; exit 1; }
done
rm -rf "$W"
echo "DONE $tag"

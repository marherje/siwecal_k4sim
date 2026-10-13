#!/bin/bash
# Prepare (and, with SUBMIT=1, submit) analysis/run_sim_chain_final.sh on condor: one job per group of
# analysis/final_v6_groups.txt whose ddsim chunks are all on EOS and whose final_v6 output (digi + nodigi
# valtrees) is incomplete.
#
#   bash analysis/submit_sim_chain_final.sh [<dir for the .sub/logs>]
#   SUBMIT=1 bash analysis/submit_sim_chain_final.sh
#
# NCHUNKS: chunks a complete group has (default: NEVENTS / CHUNK_EVENTS of launch_beam_final.sh = 50000 / 2000 = 25).
# GROUPS_FILE, PARTICLE (e-), IN_DIR (where the ddsim chunks are) as for launch_beam_final.sh.
# ALL=1 FORCE_DIGI=1 : every group, digitised branch redone (after a digitiser change). ONLY_SET=th230 restricts.
R="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
D=${1:-$R/analysis/condor_sim_chain_final}
G=${IN_DIR:-/eos/experiment/drdcalo/siw-ecal/TB2026-06/Simulation/Generated/final_v6}
GROUPS_FILE=${GROUPS_FILE:-$R/analysis/final_v6_groups.txt}
P=${PARTICLE:-e-}
O=/eos/experiment/drdcalo/siw-ecal/TB2026-06/Simulation/Processed/adc_vs_tb/final_v6
NCHUNKS=${NCHUNKS:-$(( (${NEVENTS:-50000} + ${CHUNK_EVENTS:-2000} - 1) / ${CHUNK_EVENTS:-2000} ))}
mkdir -p "$D/logs"
list="$D/jobs_$(date +%Y%m%d_%H%M%S).txt"; : > "$list"
nwait=0
while read -r th pos E xy runs; do
  [ -n "${ONLY_SET:-}" ] && [ "$th" != "$ONLY_SET" ] && continue
  case "$P" in e-) tag=${th}_${pos}_e${E} ;; *) tag=${th}_${pos}_${P%[-+]}${E} ;; esac
  n=$(ls "$G"/output_beam_${P}_${E}GeV_${pos}_xy_*_c[0-9][0-9][0-9].edm4hep.root 2>/dev/null | wc -l)
  if [ "$n" -lt "$NCHUNKS" ]; then echo "$tag: $n/$NCHUNKS ddsim chunks, waiting"; nwait=$((nwait + 1)); continue; fi
  done_=1
  for k in digi nodigi; do ls "$O/pid/${tag}_$k/"*.valtree.root >/dev/null 2>&1 || done_=0; done
  [ $done_ = 1 ] && [ "${ALL:-0}" != 1 ] || echo "$th $pos $E" >> "$list"
done < <(grep -v '^#' "$GROUPS_FILE" | awk 'NF>=5')
n=$(wc -l < "$list")
echo "--- $n groups to process, $nwait waiting for ddsim -> $list"
[ "$n" -gt 0 ] || { rm -f "$list"; exit 0; }
cat > "$D/sim_chain_final.sub" <<EOF
universe                = vanilla
executable              = /bin/bash
arguments               = $R/analysis/run_sim_chain_final.sh \$(TH) \$(POS) \$(E)
environment             = "FORCE_DIGI=${FORCE_DIGI:-0} FORCE=${FORCE:-0} PARTICLE=$P IN_DIR=$G"
log                     = $D/logs/sim_chain_final.log
output                  = $D/logs/sim_chain_\$(TH)_\$(POS)_e\$(E).out
error                   = $D/logs/sim_chain_\$(TH)_\$(POS)_e\$(E).err
request_cpus            = 1
request_memory          = 8000
request_disk            = 20000M
should_transfer_files   = NO
getenv                  = False
+JobFlavour             = "tomorrow"
queue TH, POS, E from $list
EOF
echo "--- submit file: $D/sim_chain_final.sub"
if [ "${SUBMIT:-0}" = 1 ]; then
  condor_submit "$D/sim_chain_final.sub"
else
  echo "--- not submitted (SUBMIT=1 to submit, or: condor_submit $D/sim_chain_final.sub)"
fi

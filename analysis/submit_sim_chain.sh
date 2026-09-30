#!/bin/bash
# Submit analysis/run_sim_chain.sh to condor, one job per group of analysis/final_v5_groups.txt whose ddsim sample
# exists and whose final_v5 output (digi + nodigi PID and valtrees) is incomplete.
#
#   bash analysis/submit_sim_chain.sh [<out_dir for the .sub/logs>]
#
# Re-run it after more ddsim samples land: groups already complete are skipped.
# ALL=1 FORCE_DIGI=1 bash analysis/submit_sim_chain.sh : every group, digitised branch redone (after a digitiser change).
# ONLY_SET=th230 restricts to one threshold set.
R="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
D=${1:-$R/analysis/condor_sim_chain}
G=/eos/experiment/drdcalo/siw-ecal/TB2026-06/Simulation/Generated
O=/eos/experiment/drdcalo/siw-ecal/TB2026-06/Simulation/Processed/adc_vs_tb/final_v5
mkdir -p "$D/logs"
list="$D/jobs_$(date +%Y%m%d_%H%M%S).txt"; : > "$list"
grep -v '^#' "$R/analysis/final_v5_groups.txt" | awk 'NF>=4{print $1, $2, $3}' | while read -r th E pos; do
  x=${pos%_*}; y=${pos#*_}; tag=${th}_e${E}_${pos}
  f=$G/output_beam_e-_${E}GeV_xy_${x}_${y}_sigx13.75_sigy8.25_sigE0.02_beamline_real_w1924_mat_emz.edm4hep.root
  [ -s "$f" ] || continue
  done_=1
  for k in digi nodigi; do
    ls "$O/pid/${tag}_$k/"*.valtree.root >/dev/null 2>&1 || done_=0
  done
  [ -n "${ONLY_SET:-}" ] && [ "$th" != "$ONLY_SET" ] && continue
  [ $done_ = 1 ] && [ "${ALL:-0}" != 1 ] || echo "$E $x $y $th" >> "$list"
done
n=$(wc -l < "$list")
[ "$n" -gt 0 ] || { echo "nothing to submit"; exit 0; }
cat > "$D/sim_chain.sub" <<EOF
universe                = vanilla
executable              = /bin/bash
arguments               = $R/analysis/run_sim_chain.sh \$(E) \$(X) \$(Y) \$(TH)
log                     = $D/logs/sim_chain.log
output                  = $D/logs/sim_chain_\$(TH)_e\$(E)_\$(X)_\$(Y).out
error                   = $D/logs/sim_chain_\$(TH)_e\$(E)_\$(X)_\$(Y).err
request_cpus            = 1
request_memory          = 4000
request_disk            = 8000M
should_transfer_files   = NO
getenv                  = False
environment             = "FORCE_DIGI=${FORCE_DIGI:-0}"
max_retries             = 2
+JobFlavour             = "longlunch"
queue E, X, Y, TH from $list
EOF
condor_submit "$D/sim_chain.sub"
echo "$n job(s): $list"

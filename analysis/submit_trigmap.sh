#!/bin/bash
# Trigger-map study (analysis/trigger_maps.py): the same ddsim sample of th210 S1 at 20 and 74 GeV through
# analysis/run_sim_chain_final.sh twice, with the per-slab trigger model of final_v6 (TRIGGER_TABLES=0) and with the
# per-chip/per-channel muon trigger tables at their own absolute scale (TRIGGER_TABLES=1). Same digitiser seed in both.
#
#   SUBMIT=1 bash analysis/submit_trigmap.sh
#
# IN_DIR   ddsim chunks (default /eos/user/m/marquezh/siwecal_sim/Generated/trigmap: launch_beam_final.sh with
#          NEVENTS=10000 ONLY_SET=th210 ONLY_POS=S1 ONLY_E=20|74 OUT_DIR=<that>)
# Output   Simulation/Processed/adc_vs_tb/trigmap_{perslab,muon}/{trees,pid}/
R="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
D=$R/analysis/condor_trigmap
G=${IN_DIR:-/eos/user/m/marquezh/siwecal_sim/Generated/trigmap}
O=/eos/experiment/drdcalo/siw-ecal/TB2026-06/Simulation/Processed/adc_vs_tb
mkdir -p "$D/logs"
list=$D/jobs.txt; : > "$list"
for E in 20 74; do
  n=$(ls "$G"/output_beam_e-_${E}GeV_S1_xy_*_c[0-9][0-9][0-9].edm4hep.root 2>/dev/null | wc -l)
  [ "$n" -ge 5 ] || { echo "e$E: $n/5 ddsim chunks, not ready" >&2; exit 1; }
  echo "$E 0 $O/trigmap_perslab" >> "$list"
  echo "$E 1 $O/trigmap_muon" >> "$list"
done
cat > "$D/trigmap.sub" <<EOF
universe                = vanilla
executable              = /bin/bash
arguments               = $R/analysis/run_sim_chain_final.sh th210 S1 \$(E) \$(OUT)
environment             = "TRIGGER_TABLES=\$(TT) PARTICLE=e- IN_DIR=$G"
log                     = $D/logs/trigmap.log
output                  = $D/logs/trigmap_e\$(E)_tt\$(TT).out
error                   = $D/logs/trigmap_e\$(E)_tt\$(TT).err
request_cpus            = 1
request_memory          = 8000
request_disk            = 20000M
should_transfer_files   = NO
getenv                  = False
+JobFlavour             = "tomorrow"
queue E, TT, OUT from $list
EOF
cat "$list"
if [ "${SUBMIT:-0}" = 1 ]; then condor_submit "$D/trigmap.sub"; else echo "--- not submitted: condor_submit $D/trigmap.sub"; fi

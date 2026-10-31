#!/bin/bash
# Batch job: every plot of the final_v9 analysis (border cut 50 mm), into <out>/{nocut,b50,control_v7,validation}.
#   nocut       final_panels.py without the border cut: the MIP -> GeV line (mip_gev_line.txt) and the reference panels
#   b50         final_panels.py --border-cut 50 (panels per group, resolution/linearity, border_survival.txt), the
#               summarized panels (final_panels_select.py) and the longitudinal profiles (final_profiles.py)
#   control_v7  final_v7 against the current Reconstruction without the cut: must reproduce the final-sample report
#   validation  validate_v8.py for final_v9 (hits, E/hit, profile chi2, R90 against data and sim_v7)
# Data: Reconstruction, whose PID outputs carry barxy_min_border_distance since 2026-10-31
# (replace_reconstruction_pid.sh); DATA overrides it.
#
#   bash analysis/final_v9_plots_job.sh <out>
OUT=${1:?out}
R="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
B=/eos/experiment/drdcalo/siw-ecal/TB2026-06
DATA=${DATA:-$B/Reconstruction}
SIM=$B/Simulation/Processed/final_v9
source "$R/init_key4hep.sh" >/dev/null 2>&1
export MPLCONFIGDIR=${TMPDIR:-/tmp}/mpl_$USER
cd "$R" || exit 1
for d in "$OUT" "$OUT/nocut" "$OUT/b50" "$OUT/control_v7" "$OUT/validation"; do [ -d "$d" ] || mkdir "$d"; done
step() { echo "=== $*"; "$@" || { echo "FAILED: $*"; exit 1; }; }

export FINAL_SIM=$SIM FINAL_RECO=$DATA FINAL_RECO_TREES=$B/Reconstruction
step python3 analysis/final_panels.py "$OUT/nocut"
CAL=$OUT/nocut/mip_gev_line.txt
step python3 analysis/final_panels.py "$OUT/b50" --border-cut 50 --cal-file "$CAL"
for g in th210_e20_26_-52 th210_e74_22_-49 th230_E4_e20 th230_E4_e74; do
  step python3 analysis/final_panels_select.py "$OUT/b50" --group "$g" --border-cut 50 --cal-file "$CAL"
done
step python3 analysis/final_profiles.py "$OUT/b50" --groups th210_e20_26_-52,th210_e74_22_-49,th230_E4_e20,th230_E4_e74 \
  --border-cut 50 --cal-file "$CAL"

FINAL_SIM=$B/Simulation/Processed/final_v7 FINAL_RECO=$B/Reconstruction \
  step python3 analysis/final_panels.py "$OUT/control_v7"

step python3 analysis/validate_v8.py "$OUT/validation" --new "$SIM" --cache-tag final_v9 \
  --label "sim_v9 (H2-like line + measured sag map)"
echo "DONE"

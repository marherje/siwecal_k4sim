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
for d in "$OUT" "$OUT/nocut" "$OUT/control_v7" "$OUT/validation"; do [ -d "$d" ] || mkdir "$d"; done
step() { echo "=== $*"; "$@" || { echo "FAILED: $*"; exit 1; }; }

# MIP -> GeV: one straight line per threshold set where final_panels.SET_LINES gives a scan (th230 on E2), the th210
# scan line for the others. CUTS: border cuts to draw [mm] (default 50); CUTS_ONLY=1 skips the no-cut pass, the sim_v7 control and the
# validation and reads the MIP -> GeV line of an earlier no-cut pass in <out>/nocut.
CUTS=${CUTS:-50}
export FINAL_SIM=$SIM FINAL_RECO=$DATA FINAL_RECO_TREES=$B/Reconstruction
[ "${CUTS_ONLY:-0}" = 1 ] || step python3 analysis/final_panels.py "$OUT/nocut"
CAL=$OUT/nocut/mip_gev_line.txt
for c in $CUTS; do
  D=$OUT/b$c; [ -d "$D" ] || mkdir "$D"
  step python3 analysis/final_panels.py "$D" --border-cut "$c" --cal-file "$CAL"
  for g in th210_e20_26_-52 th210_e74_22_-49 th230_E4_e20 th230_E4_e74; do
    step python3 analysis/final_panels_select.py "$D" --group "$g" --border-cut "$c" --cal-file "$CAL"
  done
  step python3 analysis/final_profiles.py "$D" --groups th210_e20_26_-52,th210_e74_22_-49,th230_E4_e20,th230_E4_e74 \
    --border-cut "$c" --cal-file "$CAL"
  step python3 analysis/contained_resolution.py "$D" "$OUT/nocut" --cut "$c" --margin 5
  step python3 analysis/position_scans.py "$D" --title "shower-like events, barxy_min_border_distance > $c mm, final_v9"
done
# SKIP_CONTROL=1: no sim_v7 control pass and no validation (neither depends on the calibration of the panels)
[ "${CUTS_ONLY:-0}" = 1 ] || [ "${SKIP_CONTROL:-0}" = 1 ] && { echo "DONE"; exit 0; }

FINAL_SIM=$B/Simulation/Processed/final_v7 FINAL_RECO=$B/Reconstruction \
  step python3 analysis/final_panels.py "$OUT/control_v7"

step python3 analysis/validate_v8.py "$OUT/validation" --new "$SIM" --cache-tag final_v9 \
  --label "sim_v9 (H2-like line + measured sag map)"
echo "DONE"

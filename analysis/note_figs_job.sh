#!/bin/bash
# Batch job: redraw every figure of the upstream-line / chip-occupancy note (beamline_occupancy_figs.py,
# validate_v8.py) into <figdir>, regenerating the 74 GeV test samples first (note_test_samples.sh) and keeping their
# ecal trees in <sampledir> so later redraws need no new simulation.
#
#   bash analysis/note_figs_job.sh <figdir> <sampledir>
FIG=${1:?figdir}; KEEP=${2:?sampledir}
R="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
W=${TMPDIR:-/tmp}/note_figs_$USER
mkdir -p "$W"
bash "$R/analysis/note_test_samples.sh" "$W"
for f in t2/d_0 t2/d_on t2/d_nooff t3/d_0; do
  [ -s "$W/$f/ecal.root" ] || { echo "missing $W/$f/ecal.root"; tail -20 "$W/$f/job3.log"; exit 1; }
done
for d in "$KEEP" "$KEEP/t2" "$KEEP/t3" "$KEEP/t2/d_0" "$KEEP/t2/d_on" "$KEEP/t2/d_nooff" "$KEEP/t3/d_0" "$FIG"; do
  [ -d "$d" ] || mkdir "$d"
done
for f in t2/d_0 t2/d_on t2/d_nooff t3/d_0; do cp "$W/$f/ecal.root" "$W/$f/job3.log" "$KEEP/$f/"; done
source "$R/init_key4hep.sh" >/dev/null 2>&1
export MPLCONFIGDIR=$W/mpl
cd "$R" || exit 1
python3 analysis/beamline_occupancy_figs.py "$FIG" --test-dir "$W" || { echo "beamline_occupancy_figs FAILED"; exit 1; }
python3 analysis/validate_v8.py "$FIG" || { echo "validate_v8 FAILED"; exit 1; }
ls -la "$FIG"
echo "DONE"

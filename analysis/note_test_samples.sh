#!/bin/bash
# Test samples of the upstream-line / chip-occupancy note (delay.png, impl_check.png in beamline_occupancy_figs.py):
# 74 GeV e- at (0, -32) mm, th230.
#   t2: H2-like line (SND_compact_beamline.xml, beam at -65.5 m), 200 events, digitised three times:
#       d_0 (OccupancySlope 0), d_on (default k*Q), d_nooff (BEAM_Z_MM=-2650: no TimeOffsetNs)
#   t3: final_v7 geometry (SND_compact_beamline_v7.xml, beam at -2.65 m), 100 events, d_0
#
#   bash analysis/note_test_samples.sh <outdir>
# No set -e / set -u: key4hep's setup and k4run do not survive them.
OUT=${1:?outdir}
R="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

digitise() {   # <sample dir> <tag> [extra env]
  local d=$1/d_$2; shift 2
  mkdir -p "$d"
  (cd "$d" && env "$@" CALIB_THRESHOLD=th230 DIGI_MODE=real RECO_TABLE_THRESHOLD=th210 INPUT_FILE=../sim.edm4hep.root \
    ECAL_TREE_OUTPUT=ecal.root k4run "$R/gaudi_jobs/pid2026_common/job3_digitize.py" > job3.log 2>&1)
  echo "$d rc=$?"; grep "beam origin\|wrote" "$d/job3.log"
}

mkdir -p "$OUT/t2" "$OUT/t3"
(cd "$OUT/t2" && COMPACT_FILE=SND_compact_beamline.xml BEAM_Z_MM=-65500 PHYSICS_LIST=QGSP_BERT_EMZ KEEP_LOCAL=1 \
  bash "$R/simulation/run_script/ddsim_beam_chunk.sh" "$R" 200 e- 74 0 -32 10.7 5.9 0.02 777 - > ddsim.log 2>&1)
(cd "$OUT/t3" && COMPACT_FILE=SND_compact_beamline_v7.xml BEAM_Z_MM=-2650 PHYSICS_LIST=QGSP_BERT_EMZ KEEP_LOCAL=1 \
  bash "$R/simulation/run_script/ddsim_beam_chunk.sh" "$R" 100 e- 74 0 -32 10.7 5.9 0.02 778 - > ddsim.log 2>&1)
source "$R/init_key4hep.sh" >/dev/null 2>&1
export LD_LIBRARY_PATH=$R/install/lib64:$R/install/lib:$LD_LIBRARY_PATH
export PYTHONPATH=$R/install/lib64:$R/install/lib:$R/install/python:$PYTHONPATH
digitise "$OUT/t2" 0 OCCUPANCY_SLOPE=0 BEAM_Z_MM=-65500
digitise "$OUT/t2" on BEAM_Z_MM=-65500
digitise "$OUT/t2" nooff BEAM_Z_MM=-2650
digitise "$OUT/t3" 0 OCCUPANCY_SLOPE=0 BEAM_Z_MM=-2650

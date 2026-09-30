#!/bin/bash
# One beam sample through the analysis chain: ddsim edm4hep -> digitised (real, reconstructed with the th210 table)
# and undigitised (simple) ecal trees -> ECalPid edm4hep (hit MIP cut 0.5).
#
#   bash analysis/run_sim_chain.sh <E_GeV> <x_mm> <y_mm> <set: th210|th220|th230> [<out_base>]
#   (also writes the valtree next to each PID file; condor: analysis/submit_sim_chain.sh)
#
# Input : Simulation/Generated/output_beam_e-_<E>GeV_xy_<x>_<y>_sigx13.75_sigy8.25_sigE0.02_beamline_real_w1924_mat_emz.edm4hep.root
# Output: <out_base>/trees/ecal_<set>_e<E>_<x>_<y>{,_simple}.root
#         <out_base>/pid/<set>_e<E>_<x>_<y>_{digi,nodigi}/*.edm4hep.root
# out_base defaults to Simulation/Processed/adc_vs_tb/final_v5.
#
# Environment:
#   DIGI_CALIB_DIR   MuonCalib_gaudi-style table tree for the digitiser (default: masking_info/calibration/MuonCalib_gaudi_fixed)
#   RECO_TABLE       table that reconstructs the simulated ADC (default th210, like Reconstructed_final)
#   TB_REPO          siwecal-tb2026 checkout for the PID stage (its gaudi_source/build must carry fractal_dimension)
#   WORK             scratch dir (default ${TMPDIR:-/tmp}/sim_chain_$USER)
#   PED_NOISE        1 (default): per-cell noise and true gain spread in the digitiser; 0 = the scalar noise
#   FORCE=1          redo steps whose output exists
# Rewritten from final_v3's chain_pos.sh + digi_recoth210.sh.
E=${1:?energy}; X=${2:?x}; Y=${3:?y}; SET=${4:?set}
R="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
T=${TB_REPO:-$(cd "$R/../siwecal-tb2026" && pwd)}
B=/eos/experiment/drdcalo/siw-ecal/TB2026-06/Simulation
O=${5:-$B/Processed/adc_vs_tb/final_v5}
export DIGI_CALIB_DIR=${DIGI_CALIB_DIR:-$R/masking_info/calibration/MuonCalib_gaudi_fixed}
RECO_TABLE=${RECO_TABLE:-th210}
# per-cell electronic noise from the fixed pedestal tables + gain spread without fit error (job3_digitize.py)
export PED_NOISE=${PED_NOISE:-1}
in=$B/Generated/output_beam_e-_${E}GeV_xy_${X}_${Y}_sigx13.75_sigy8.25_sigE0.02_beamline_real_w1924_mat_emz.edm4hep.root
tag=${SET}_e${E}_${X}_${Y}
W=${WORK:-${TMPDIR:-/tmp}/sim_chain_$USER}/$tag

# key4hep before any `set -u` (it would leave k4run off the PATH), then the repo's own components.
source "$R/init_key4hep.sh" >/dev/null 2>&1
command -v k4run >/dev/null || { echo "k4run not on PATH" >&2; exit 1; }
set -o pipefail
export LD_LIBRARY_PATH=$R/install/lib64:$R/install/lib:$LD_LIBRARY_PATH
export PYTHONPATH=$R/install/lib64:$R/install/lib:$R/install/python:$PYTHONPATH

[ -s "$in" ] || { echo "missing ddsim input $in" >&2; exit 1; }
[ -d "$DIGI_CALIB_DIR/mips/$RECO_TABLE" ] || { echo "no $RECO_TABLE table in $DIGI_CALIB_DIR" >&2; exit 1; }
# plain mkdir per level: mkdir -p walks up into EOS parents it cannot touch
for d in "$O" "$O/trees" "$O/pid" "$O/pid/${tag}_digi" "$O/pid/${tag}_nodigi"; do [ -d "$d" ] || mkdir "$d"; done
mkdir -p "$W"

for mode in real simple; do
  if [ $mode = simple ]; then
    out=$O/trees/ecal_${tag}_simple.root; pid=$O/pid/${tag}_nodigi; reco=""
    coll=(-c SiPadHitsMapped --masking-collection SiPadHitsMasked)
  else
    out=$O/trees/ecal_${tag}.root; pid=$O/pid/${tag}_digi; reco=$RECO_TABLE
    coll=(-c SiPadHitsRealAdc --masking-collection "")
  fi
  if [ "${FORCE:-0}" = 1 ] || [ ! -s "$out" ]; then
    (cd "$W" && CALIB_THRESHOLD=$SET DIGI_MODE=$mode RECO_TABLE_THRESHOLD=$reco INPUT_FILE=$in \
      k4run "$R/gaudi_jobs/pid2026_common/job3_digitize.py" > "$W/job3_$mode.log" 2>&1) \
      || { echo "job3 FAILED $tag $mode"; tail -5 "$W/job3_$mode.log"; exit 1; }
    python3 "$R/analysis/sim_to_ecal_tree.py" -i "$W/digitized.edm4hep.root" -o "$W/ecal.root" "${coll[@]}" \
      > "$W/conv_$mode.log" 2>&1 || { echo "conv FAILED $tag $mode"; tail -5 "$W/conv_$mode.log"; exit 1; }
    cp "$W/ecal.root" "$out" && rm -f "$W/digitized.edm4hep.root" "$W/ecal.root"
  fi
  if [ "${FORCE:-0}" = 1 ] || ! ls "$pid"/*.edm4hep.root >/dev/null 2>&1; then
    (cd "$T" && source "$T/setup.sh" >/dev/null 2>&1 && \
      python3 gaudi_jobs/run_pid_batch.py --file "$out" --outdir "$W/pid_$mode" --format edm4hep --hit-mip-cut 0.5 \
      > "$W/pid_$mode.log" 2>&1) || { echo "pid FAILED $tag $mode"; tail -5 "$W/pid_$mode.log"; exit 1; }
    cp "$W/pid_$mode"/*.edm4hep.root "$pid"/
  fi
  # valtree next to the PID file (gaudi_jobs/pid_to_valtree.py skips an existing one)
  (cd "$T" && source "$T/setup.sh" >/dev/null 2>&1 && \
    python3 gaudi_jobs/pid_to_valtree.py "$pid"/*.edm4hep.root > "$W/valtree_$mode.log" 2>&1) \
    || { echo "valtree FAILED $tag $mode"; tail -5 "$W/valtree_$mode.log"; exit 1; }
done
rm -rf "$W"
echo "DONE $tag"

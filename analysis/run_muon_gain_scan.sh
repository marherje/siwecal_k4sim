#!/bin/bash
# Digitise one muon ddsim sample at a list of ADC gains for one threshold set, for analysis/muon_gain_scan.py.
#
#   bash analysis/run_muon_gain_scan.sh <muon ddsim edm4hep> <set: th210|th220|th230> <outdir> <gain> [<gain> ...]
#
# Each gain -> <outdir>/<set>_mu_g<gain>.root (ecal tree of SiPadHitsRealAdc). Same digitiser configuration as the
# final_v5 samples: DIGI_CALIB_DIR (default masking_info/calibration/MuonCalib_gaudi_fixed), per-cell pedestal
# noise (PED_NOISE=1), reconstruction with the th210 table, only ADC_PER_MIP_OVERRIDE changing.
IN=${1:?muon ddsim file}; SET=${2:?set}; OUT=${3:?outdir}; shift 3
R="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
source "$R/init_key4hep.sh" >/dev/null 2>&1
command -v k4run >/dev/null || { echo "k4run not on PATH" >&2; exit 1; }
set -o pipefail
export LD_LIBRARY_PATH=$R/install/lib64:$R/install/lib:$LD_LIBRARY_PATH
export PYTHONPATH=$R/install/lib64:$R/install/lib:$R/install/python:$PYTHONPATH
export DIGI_CALIB_DIR=${DIGI_CALIB_DIR:-$R/masking_info/calibration/MuonCalib_gaudi_fixed}
export PED_NOISE=${PED_NOISE:-1}
[ -d "$OUT" ] || mkdir "$OUT"
for g in "$@"; do
  out=$OUT/${SET}_mu_g$g.root
  [ -s "$out" ] && [ "${FORCE:-0}" != 1 ] && { echo "exists $out"; continue; }
  W=$(mktemp -d "${TMPDIR:-/tmp}/gainscan_XXXXXX")
  (cd "$W" && CALIB_THRESHOLD=$SET DIGI_MODE=real RECO_TABLE_THRESHOLD=th210 ADC_PER_MIP_OVERRIDE=$g INPUT_FILE=$IN \
     k4run "$R/gaudi_jobs/pid2026_common/job3_digitize.py" > "$W/job3.log" 2>&1) \
    || { echo "job3 FAILED g=$g"; tail -5 "$W/job3.log"; exit 1; }
  (cd "$R" && python3 -m analysis.sim_to_ecal_tree -i "$W/digitized.edm4hep.root" -o "$W/ecal.root" \
     -c SiPadHitsRealAdc --masking-collection SiPadHitsRealMasked > "$W/conv.log" 2>&1) \
    || { echo "conv FAILED g=$g"; tail -5 "$W/conv.log"; exit 1; }
  cp "$W/ecal.root" "$out" && rm -rf "$W"
  echo "DONE $out"
done

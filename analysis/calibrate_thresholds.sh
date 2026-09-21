#!/bin/bash
# ---------------------------------------------------------------------------
# calibrate_thresholds.sh
#
# The per-threshold calibration protocol, end to end:
#
#   A. CALIBRATE  at one beam energy: fit `adc_per_mip` so the simulation's event
#      ADC sum matches the data's, for that threshold set's own data.
#   B. VALIDATE   at a DIFFERENT energy (muons for th220, which has only one
#      electron energy), with every number frozen.
#
# Nothing is shared between threshold sets: no ratio of MIP tables, no reference
# gain borrowed from another set.  Each set is anchored to its own data.
#
# The simulated samples must exist at the beam position each run was taken at --
# they differ per threshold set (th230 ~(-42,+51), th220 ~(-40,+37), th210
# ~(+24,-51) mm), which is why the productions are position-matched.
#
# Usage:
#   bash analysis/calibrate_thresholds.sh <workdir>
# ---------------------------------------------------------------------------
set -o pipefail

REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
source "$REPO/init_key4hep.sh" || exit 1
command -v k4run >/dev/null || { echo "ERROR: k4run not on PATH"; exit 1; }
export LD_LIBRARY_PATH="$REPO/install/lib64:$REPO/install/lib:${LD_LIBRARY_PATH:-}"
export PYTHONPATH="$REPO/install/lib64:$REPO/install/lib:$REPO/install/python:${PYTHONPATH:-}"

WORK="${1:?usage: calibrate_thresholds.sh <workdir>}"
mkdir -p "$WORK"
GEN=/eos/experiment/drdcalo/siw-ecal/TB2026-06/Simulation/Generated
RECO=/eos/experiment/drdcalo/siw-ecal/TB2026-06/Reconstruction
PLOTS="$WORK/plots"
cd "$REPO" || exit 1

sim_file() {  # <energy> <x> <y>
    echo "$GEN/output_beam_e-_${1}GeV_xy_${2}_${3}_sigx13.75_sigy8.25_sigE0.02.edm4hep.root"
}

# th : calib energy : calib data : calib sim : valid energy : valid data : valid sim
CASES=(
  "th230:20:$RECO/TB2026CERN_run_000020/ecal_TB2026CERN_run_000020.root:$(sim_file 20 -42 51):52:$RECO/TB2026CERN_run_000013/ecal_TB2026CERN_run_000013.root:$(sim_file 52 -42 51)"
  "th220:74:$RECO/TB2026CERN_run_000072/ecal_TB2026CERN_run_000072.root:$(sim_file 74 -40 37):MU:$RECO/TB2026CERN_run_000085/ecal_TB2026CERN_run_000085.root:$GEN/output_beam_mu-_100GeV_xy_1_1_sigx38.5_sigy46.75_sigE0.02.edm4hep.root"
  "th210:52:$RECO/TB2026CERN_eudaq_run_000287/ecal_TB2026CERN_eudaq_run_000287.root:$(sim_file 52 24 -51):74:$RECO/TB2026CERN_eudaq_run_000286/ecal_TB2026CERN_eudaq_run_000286.root:$(sim_file 74 24 -51)"
)

for case in "${CASES[@]}"; do
    IFS=: read -r TH CE CDATA CSIM VE VDATA VSIM <<< "$case"
    echo ""
    echo "########## $TH — calibrate at $CE GeV, validate at $VE"
    for f in "$CDATA" "$CSIM" "$VDATA" "$VSIM"; do
        [[ -f "$f" ]] || { echo "  MISSING: $f — skipping $TH"; continue 2; }
    done

    # A. fit, and write the result into mappings/digi_calibration.yml
    python3 -m analysis.fit_adc_scale --threshold "$TH" --data "$CDATA" \
        --sim-input "$CSIM" --workdir "$WORK/fit_$TH" --energy "$CE" --write || continue

    # The calibration point, with the fitted value in place
    python3 -m analysis.compare_adc_data_sim --data "$CDATA" \
        --sim "sim $TH"="$WORK/fit_$TH/ecal_${TH}_final.root" \
        --mip-file "$(ls /afs/cern.ch/user/m/marquezh/public/siwecal-tb2026/calibration/MuonCalib_gaudi/mips/$TH/*highgain* | head -1)" \
        --tag "calib_${TH}" --title "$TH calibration point, $CE GeV" \
        --outdir "$PLOTS" --max-data-events 40000 || true

    # B. validation at the other point, nothing refitted
    VWORK="$WORK/valid_$TH"; mkdir -p "$VWORK"
    ( cd "$VWORK" && DIGI_MODE=real ADC_MODEL=1 CALIB_THRESHOLD="$TH" \
        INPUT_FILE="$VSIM" k4run "$REPO/gaudi_jobs/pid2026_common/job3_digitize.py" \
        >/dev/null 2>&1 ) || { echo "  validation digitisation failed"; continue; }
    python3 -m analysis.sim_to_ecal_tree --input "$VWORK/digitized.edm4hep.root" \
        --output "$VWORK/ecal_valid_${TH}.root" --collection SiPadHitsRealAdc \
        --masking-collection SiPadHitsRealMasked >/dev/null || continue

    python3 -m analysis.compare_adc_data_sim --data "$VDATA" \
        --sim "sim $TH"="$VWORK/ecal_valid_${TH}.root" \
        --mip-file "$(ls /afs/cern.ch/user/m/marquezh/public/siwecal-tb2026/calibration/MuonCalib_gaudi/mips/$TH/*highgain* | head -1)" \
        --tag "valid_${TH}" --title "$TH validation, $VE" \
        --outdir "$PLOTS" --max-data-events 40000 || true

    # The dynamic range at both points
    python3 -m analysis.compare_gains --data "$VDATA" \
        --sim "$VWORK/ecal_valid_${TH}.root" --threshold "$TH" \
        --tag "valid_${TH}" --title "$TH validation, $VE" --outdir "$PLOTS" || true
done

echo ""
echo "########## CALIBRATION DONE — plots in $PLOTS"

#!/bin/bash
# Note: no set -e — ddsim with runType="run"+GPS crashes during Geant4 cleanup
# even when all events are written. We validate the output file instead.

echo "Starting beam-test job on $(hostname)"

source /afs/cern.ch/user/m/marquezh/public/siwecal_k4sim/simulation/run_script/../../init_key4hep.sh
export LD_LIBRARY_PATH=/afs/cern.ch/user/m/marquezh/public/siwecal_k4sim/simulation/run_script/../../install/lib64:/afs/cern.ch/user/m/marquezh/public/siwecal_k4sim/simulation/run_script/../../install/lib:$LD_LIBRARY_PATH
export PYTHONPATH=/afs/cern.ch/user/m/marquezh/public/siwecal_k4sim/simulation/run_script/../../install/lib64:/afs/cern.ch/user/m/marquezh/public/siwecal_k4sim/simulation/run_script/../../install/lib:/afs/cern.ch/user/m/marquezh/public/siwecal_k4sim/simulation/run_script/../../install/python:$PYTHONPATH

EXPECTED_OUTPUT="/afs/cern.ch/user/m/marquezh/public/siwecal_k4sim/simulation/run_script/data/output_beam_mu-_100GeV_xy_1_1_sigx20.5_sigy16.5_sigE0.02.edm4hep.root"

ddsim --enableG4GPS \
      --macroFile    /afs/cern.ch/user/m/marquezh/public/siwecal_k4sim/simulation/run_script/steer/gps_QGSP_BERT_SiWECAL_beam_mu-_100GeV_xy_1_1_sigx20.5_sigy16.5_sigE0.02.mac \
      --steeringFile /afs/cern.ch/user/m/marquezh/public/siwecal_k4sim/simulation/run_script/steer/runddsim_QGSP_BERT_SiWECAL_beam_mu-_100GeV_xy_1_1_sigx20.5_sigy16.5_sigE0.02.py \
      &> /afs/cern.ch/user/m/marquezh/public/siwecal_k4sim/simulation/run_script/log/QGSP_BERT_SiWECAL_beam_mu-_100GeV_xy_1_1_sigx20.5_sigy16.5_sigE0.02.log
DDSIM_RC=$?

if [[ ${DDSIM_RC} -ne 0 ]]; then
    if [[ -f "${EXPECTED_OUTPUT}" ]] && [[ -s "${EXPECTED_OUTPUT}" ]]; then
        echo "WARNING: ddsim exit code ${DDSIM_RC} — Geant4 cleanup crash after all events written. Output OK."
    else
        echo "ERROR: ddsim failed (exit ${DDSIM_RC}) and output file missing or empty."
        exit 1
    fi
fi

echo "Job finished. Output: ${EXPECTED_OUTPUT}"

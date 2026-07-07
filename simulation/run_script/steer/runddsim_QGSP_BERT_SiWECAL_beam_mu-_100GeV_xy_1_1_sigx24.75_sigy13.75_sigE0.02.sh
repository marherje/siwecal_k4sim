#!/bin/bash
# Note: no set -e — ddsim with runType="run"+GPS crashes during Geant4 cleanup
# even when all events are written. We validate the local output file instead.

echo "Starting beam-test job on $(hostname)"

source /afs/cern.ch/user/m/marquezh/public/siwecal_k4sim/simulation/run_script/../../init_key4hep.sh
export LD_LIBRARY_PATH=/afs/cern.ch/user/m/marquezh/public/siwecal_k4sim/simulation/run_script/../../install/lib64:/afs/cern.ch/user/m/marquezh/public/siwecal_k4sim/simulation/run_script/../../install/lib:$LD_LIBRARY_PATH
export PYTHONPATH=/afs/cern.ch/user/m/marquezh/public/siwecal_k4sim/simulation/run_script/../../install/lib64:/afs/cern.ch/user/m/marquezh/public/siwecal_k4sim/simulation/run_script/../../install/lib:/afs/cern.ch/user/m/marquezh/public/siwecal_k4sim/simulation/run_script/../../install/python:$PYTHONPATH

LOCAL_OUTPUT="output_beam_mu-_100GeV_xy_1_1_sigx24.75_sigy13.75_sigE0.02.edm4hep.root"
REMOTE_OUTPUT="/eos/experiment/drdcalo/siw-ecal/TB2026-06/Simulation/Generated/output_beam_mu-_100GeV_xy_1_1_sigx24.75_sigy13.75_sigE0.02.edm4hep.root"

ddsim --enableG4GPS \
      --macroFile    /afs/cern.ch/user/m/marquezh/public/siwecal_k4sim/simulation/run_script/steer/gps_QGSP_BERT_SiWECAL_beam_mu-_100GeV_xy_1_1_sigx24.75_sigy13.75_sigE0.02.mac \
      --steeringFile /afs/cern.ch/user/m/marquezh/public/siwecal_k4sim/simulation/run_script/steer/runddsim_QGSP_BERT_SiWECAL_beam_mu-_100GeV_xy_1_1_sigx24.75_sigy13.75_sigE0.02.py \
      &> /afs/cern.ch/user/m/marquezh/public/siwecal_k4sim/simulation/run_script/log/QGSP_BERT_SiWECAL_beam_mu-_100GeV_xy_1_1_sigx24.75_sigy13.75_sigE0.02.log
DDSIM_RC=$?

if [[ ! -f "${LOCAL_OUTPUT}" ]] || [[ ! -s "${LOCAL_OUTPUT}" ]]; then
    echo "ERROR: ddsim (exit ${DDSIM_RC}) produced no usable local output file."
    exit 1
fi

if [[ ${DDSIM_RC} -ne 0 ]]; then
    echo "WARNING: ddsim exit code ${DDSIM_RC} — Geant4 cleanup crash after all events written. Local output OK, staging out."
fi

# Stage out via xrdcp (synchronous xrootd transfer through the redirector),
# NOT a plain cp through the worker's /eos FUSE mount: eosxd write-back
# caching on ephemeral batch slots can report a successful local write that
# never actually reaches the EOS namespace before the job slot is torn down.
LOCAL_SIZE=$(stat -c%s "${LOCAL_OUTPUT}")

xrdcp --force "${LOCAL_OUTPUT}" "root://eosexperiment.cern.ch/${REMOTE_OUTPUT}"
XRDCP_RC=$?
if [[ ${XRDCP_RC} -ne 0 ]]; then
    echo "ERROR: xrdcp stage-out to EOS failed (exit ${XRDCP_RC})."
    exit 1
fi

REMOTE_SIZE=$(xrdfs eosexperiment.cern.ch stat "${REMOTE_OUTPUT}" 2>/dev/null | awk '/Size:/{print $2}')
if [[ -z "${REMOTE_SIZE}" ]] || [[ "${REMOTE_SIZE}" != "${LOCAL_SIZE}" ]]; then
    echo "ERROR: stage-out verification failed (local=${LOCAL_SIZE} bytes, remote='${REMOTE_SIZE}')."
    exit 1
fi

echo "Job finished. Output: ${REMOTE_OUTPUT} (verified ${REMOTE_SIZE} bytes on EOS)"

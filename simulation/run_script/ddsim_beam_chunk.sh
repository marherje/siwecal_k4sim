#!/bin/bash
# One ddsim chunk of a GPS beam sample, run on the condor worker (executable of launch_beam_final.sh's .sub).
# Writes the GPS macro and the ddsim steering in the job's scratch dir, simulates, and stages the output out to
# EOS with xrdcp + size check (a plain cp through the worker's /eos FUSE mount can report success for a file
# that never reaches the namespace). Nothing downstream: digitisation and PID run per group afterwards
# (analysis/submit_sim_chain_final.sh), so a digitiser change never needs new Geant4.
#
#   ddsim_beam_chunk.sh <repo_root> <nevents> <particle> <E_GeV|Emin-Emax> <x_mm> <y_mm> <sigma_x_mm> <sigma_y_mm> \
#                       <sigma_E_frac> <seed> <remote_output>
#
# Environment (set in the .sub): COMPACT_FILE, BEAM_Z_MM, PHYSICS_LIST; KEEP_LOCAL=1 skips the stage-out.
# No set -e: ddsim with runType="run" + GPS can crash in the Geant4 cleanup after writing every event; the
# output file is validated instead.

repo=$1; nevents=$2; particle=$3; energy=$4; pos_x=$5; pos_y=$6; sigma_x=$7; sigma_y=$8; sigma_E=$9
seed=${10}; remote=${11}
: "${COMPACT_FILE:?}" "${BEAM_Z_MM:?}" "${PHYSICS_LIST:?}"

echo "ddsim chunk on $(hostname): ${particle} ${energy} GeV at (${pos_x}, ${pos_y}) mm, sigma (${sigma_x}, ${sigma_y}) mm," \
     "sigma_E ${sigma_E}, N=${nevents}, seed=${seed}, ${PHYSICS_LIST}, ${COMPACT_FILE}, beam z ${BEAM_Z_MM} mm"
echo "-> ${remote}"

source "${repo}/init_key4hep.sh"
command -v ddsim >/dev/null || { echo "ERROR: ddsim not on PATH"; exit 1; }
# the SiPadDetector DD4hep plugin lives in the repo's install/
export LD_LIBRARY_PATH=${repo}/install/lib64:${repo}/install/lib:${LD_LIBRARY_PATH}
export PYTHONPATH=${repo}/install/lib64:${repo}/install/lib:${repo}/install/python:${PYTHONPATH}

# <E_GeV> = "<min>-<max>" (e.g. 5-100): energy flat in [min, max] GeV (background samples); sigma_E is then ignored
if [[ "${energy}" == *-* ]]; then
  ene_cmds="/gps/ene/type Lin
/gps/ene/gradient 0
/gps/ene/intercept 1
/gps/ene/min ${energy%-*} GeV
/gps/ene/max ${energy#*-} GeV"
else
  sigma_E_GeV=$(python3 -c "print(${energy} * ${sigma_E})")
  ene_cmds="/gps/ene/type Gauss
/gps/ene/mono ${energy} GeV
/gps/ene/sigma ${sigma_E_GeV} GeV"
fi
{
echo "/gps/verbose 0"
echo "/gps/particle ${particle}"
echo "/gps/direction 0 0 1"
echo "/gps/ang/rot1 1 0 0"
echo "/gps/ang/rot2 0 1 0"
echo "/gps/pos/type Beam"
echo "/gps/pos/shape Circle"
echo "/gps/pos/centre ${pos_x} ${pos_y} ${BEAM_Z_MM} mm"
echo "/gps/pos/sigma_x ${sigma_x} mm"
echo "/gps/pos/sigma_y ${sigma_y} mm"
echo "${ene_cmds}"
echo "/run/beamOn ${nevents}"
} > gps.mac

local_out=sim.edm4hep.root
cat > steer.py <<EOF
import os
from DDSim.DD4hepSimulation import DD4hepSimulation

# DD4hep 1.35 + Python 3.13: addParametersToRunHeader returns a dict cppyy cannot convert.
try:
    from DDSim.Helper.Meta import Meta as _M
    _M.addParametersToRunHeader = lambda self, dds: {}
except Exception:
    pass

compact = os.path.abspath("${repo}/simulation/geometry/${COMPACT_FILE}")
if not os.path.isfile(compact):
    raise RuntimeError("Compact file not found: " + compact)
SIM = DD4hepSimulation()
SIM.runType = "run"          # the GPS macro drives the loop (/run/beamOn)
SIM.skipNEvents = 0
SIM.compactFile = compact
SIM._compactFile = compact
SIM.outputFile = "${local_out}"
SIM.physicsList = "${PHYSICS_LIST}"
# Per-chunk seed: without it every chunk of a point replays the same events.
SIM.random.seed = ${seed}
SIM.random.enableEventSeed = True
EOF

ddsim --enableG4GPS --macroFile gps.mac --steeringFile steer.py > ddsim.log 2>&1
rc=$?
tail -20 ddsim.log
[ -s "${local_out}" ] || { echo "ERROR: ddsim (exit ${rc}) wrote no output"; exit 1; }
n=$(python3 -c "import uproot; print(uproot.open('${local_out}')['events'].num_entries)" 2>/dev/null)
[ "${n}" = "${nevents}" ] || { echo "ERROR: ${n:-no} events in the output, expected ${nevents}"; exit 1; }
[ ${rc} -ne 0 ] && echo "WARNING: ddsim exit ${rc} after writing all ${n} events (Geant4 cleanup crash)"

# KEEP_LOCAL=1: leave sim.edm4hep.root in the job's directory for the caller (sim_digi_chunk.sh), no stage-out
if [ "${KEEP_LOCAL:-0}" = 1 ]; then echo "DONE local ${local_out} (${n} events)"; exit 0; fi
size=$(stat -c%s "${local_out}")
case "${remote}" in /eos/user/*) mgm=eosuser.cern.ch ;; *) mgm=eosexperiment.cern.ch ;; esac
xrdcp --force "${local_out}" "root://${mgm}/${remote}" || { echo "ERROR: xrdcp failed"; exit 1; }
rsize=$(xrdfs ${mgm} stat "${remote}" 2>/dev/null | awk '/Size:/{print $2}')
[ "${rsize}" = "${size}" ] || { echo "ERROR: stage-out check failed (local ${size}, remote '${rsize}')"; exit 1; }
echo "DONE ${remote} (${n} events, ${size} bytes)"

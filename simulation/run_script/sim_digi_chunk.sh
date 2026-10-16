#!/bin/bash
# One chunk of a beam point, ddsim AND digitisation in the same condor job (node of analysis/make_sim_dag.py):
# ddsim_beam_chunk.sh with KEEP_LOCAL=1 -> job3_digitize.py twice (real: digitised, reconstructed with the
# RECO_TABLE table; simple: undigitised) -> the two ecal trees staged out to EOS with xrdcp + size check.
# The raw ddsim output never leaves the worker: the whole campaign's Geant4 output (~2.7 TB at 50k events per
# point) does not fit in the project quota. Re-digitising therefore needs ddsim again.
#
#   sim_digi_chunk.sh <repo_root> <nevents> <particle> <E_GeV> <x_mm> <y_mm> <sigma_x_mm> <sigma_y_mm> \
#                     <sigma_E_frac> <seed> <set> <remote_tree_without_.root>
#   -> <remote>.root (real) and <remote>_simple.root
#
# Environment (set in the .sub): COMPACT_FILE, BEAM_Z_MM, PHYSICS_LIST, RECO_TABLE (th210); the digitiser
# defaults of job3_digitize.py (DIGI_CALIB_DIR = EOS Calib/MuonCalib, TABLE_THRESHOLD th210, PED_NOISE 1).
# No set -e / set -u: key4hep's setup and k4run do not survive them.
repo=$1; nev=$2; particle=$3; energy=$4; x=$5; y=$6; sx=$7; sy=$8; sE=$9; seed=${10}; set_=${11}; remote=${12}
: "${RECO_TABLE:=th210}"

KEEP_LOCAL=1 bash "${repo}/simulation/run_script/ddsim_beam_chunk.sh" "${repo}" "${nev}" "${particle}" "${energy}" \
  "${x}" "${y}" "${sx}" "${sy}" "${sE}" "${seed}" "-" || { echo "ERROR: ddsim chunk failed"; exit 1; }
sim=$PWD/sim.edm4hep.root

source "${repo}/init_key4hep.sh" >/dev/null 2>&1
command -v k4run >/dev/null || { echo "ERROR: k4run not on PATH"; exit 1; }
export LD_LIBRARY_PATH=${repo}/install/lib64:${repo}/install/lib:${LD_LIBRARY_PATH}
export PYTHONPATH=${repo}/install/lib64:${repo}/install/lib:${repo}/install/python:${PYTHONPATH}

case "${remote}" in /eos/user/*) mgm=eosuser.cern.ch ;; *) mgm=eosexperiment.cern.ch ;; esac
for mode in real simple; do
  if [ $mode = real ]; then reco=$RECO_TABLE; dst=${remote}.root; else reco=""; dst=${remote}_simple.root; fi
  mkdir -p "$mode"
  (cd "$mode" && CALIB_THRESHOLD=${set_} DIGI_MODE=$mode RECO_TABLE_THRESHOLD=$reco INPUT_FILE=$sim \
     ECAL_TREE_OUTPUT=ecal.root k4run "${repo}/gaudi_jobs/pid2026_common/job3_digitize.py" > job3.log 2>&1)
  rc=$?
  if [ $rc -ne 0 ]; then
    # a crash in the teardown after the tree is complete ("free(): invalid pointer") loses nothing
    if grep -q "\[EcalTreeWriter\] wrote ${nev} events" "$mode/job3.log"; then
      echo "WARNING: job3 $mode exit $rc after writing all ${nev} events"
    else
      echo "ERROR: job3 $mode failed (exit $rc)"; tail -20 "$mode/job3.log"; exit 1
    fi
  fi
  n=$(python3 -c "import uproot; print(uproot.open('$mode/ecal.root')['ecal'].num_entries)" 2>/dev/null)
  [ -n "$n" ] && [ "$n" -gt 0 ] || { echo "ERROR: empty ecal tree ($mode)"; exit 1; }
  rm -f "$mode/digitized.edm4hep.root"
  size=$(stat -c%s "$mode/ecal.root")
  xrdcp --force "$mode/ecal.root" "root://${mgm}/${dst}" || { echo "ERROR: xrdcp $mode failed"; exit 1; }
  rsize=$(xrdfs ${mgm} stat "${dst}" 2>/dev/null | awk '/Size:/{print $2}')
  [ "${rsize}" = "${size}" ] || { echo "ERROR: stage-out check failed for $dst (local ${size}, remote '${rsize}')"; exit 1; }
  echo "DONE $dst ($n events with hits, $size bytes)"
done
rm -f "$sim"

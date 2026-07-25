#!/bin/bash
# Chunked beam-test simulation + per-chunk digitisation for the SiWECAL PID
# campaign (e-/mu-/pi- at 74 and 99 GeV, 50k events per point).
#
# Difference with generic_condor_beam.sh
# --------------------------------------
#   * One condor job = one CHUNK of a point, with its own RNG seed, so a 50k
#     point is split over many jobs instead of a single ~26 h job.
#   * The job does ddsim AND the two cheap downstream steps (job3 digitisation
#     + ecal-tree conversion) on the worker's local scratch, and stages out:
#         Generated/chunks/output_<chunk_label>.edm4hep.root   (raw sim, kept
#                                                               for reprocessing)
#         Processed/chunks/<chunk_label>_ecal.root             (ecal TTree)
#     The huge digitized.edm4hep.root stays on scratch and is discarded: it can
#     always be regenerated from the raw sim chunk, and 6 points x 9 GB of it
#     would be pure dead weight.
#   * The per-point pipelines then only have to hadd the small ecal trees and
#     run k4SiWEcalReco (see gaudi_jobs/pid2026_common/pipeline_common.sh).
#
# Usage:
#   ./generic_condor_beam_chunk.sh <nevents> <particle> <energy_GeV>  \
#                                   <pos_x_mm> <pos_y_mm>             \
#                                   <sigma_x_mm> <sigma_y_mm>         \
#                                   <sigma_E_frac> <chunk_id> <seed>
#
# Set DRY_RUN=1 to generate the steer/sub files without submitting.

set -uo pipefail

mkdir -p steer log

nevents=$1
particle=$2
energy=$3
pos_x=$4
pos_y=$5
sigma_x=$6
sigma_y=$7
sigma_E=$8
chunk=$9
seed=${10}

BEAM_Z_MM=-2000
PHYSLIST="QGSP_BERT"

local=$PWD
repo_root="$(cd "${local}/../.." && pwd)"
geometry_folder="${local}/../geometry"
eos_base="/eos/experiment/drdcalo/siw-ecal/TB2026-06/Simulation"
data_path="${eos_base}/Generated/chunks"
tree_path="${eos_base}/Processed/chunks"
steer_path="${local}/steer"
log_path="${local}/log"

# Point label: identical for every chunk of the same (particle, E, position).
point=beam_${particle}_${energy}GeV_xy_${pos_x}_${pos_y}_sigx${sigma_x}_sigy${sigma_y}_sigE${sigma_E}
chunk_tag=$(printf "c%03d" "${chunk}")
label=${point}_${chunk_tag}

# Run number written to the ecal tree: keeps chunk (and energy) identity after
# the trees of several chunks - and of 74+99 GeV - are hadd'ed together.
run_number=$(( energy * 1000 + chunk ))

gpsmac=gps_${PHYSLIST}_${label}.mac
scriptname=runddsim_${PHYSLIST}_${label}.py
condorsh=runddsim_${PHYSLIST}_${label}.sh
condorsub=runddsim_${PHYSLIST}_${label}.sub
condorfile=runddsim_${PHYSLIST}_${label}

sim_name="output_${label}.edm4hep.root"
remote_sim="${data_path}/${sim_name}"
remote_tree="${tree_path}/${label}_ecal.root"

sigma_E_GeV=$(python3 -c "print(${energy} * ${sigma_E})")

echo "Chunk ${chunk_tag}: ${particle} ${energy} GeV  centre=(${pos_x},${pos_y}) mm  sigma=(${sigma_x},${sigma_y}) mm  N=${nevents}  seed=${seed}  run=${run_number}"

# --------------------------------------------------
# GEANT4 GPS MACRO FILE
# --------------------------------------------------
{
echo "/gps/verbose 0"
echo "/gps/particle ${particle}"
echo "/gps/direction 0 0 1"
echo "/gps/ang/rot1 1 0 0"
echo "/gps/ang/rot2 0 1 0"

# Spatial profile: elliptical Gaussian beam
echo "/gps/pos/type Beam"
echo "/gps/pos/shape Circle"
echo "/gps/pos/centre ${pos_x} ${pos_y} ${BEAM_Z_MM} mm"
echo "/gps/pos/sigma_x ${sigma_x} mm"
echo "/gps/pos/sigma_y ${sigma_y} mm"

# Energy distribution
if (( $(echo "${sigma_E} > 0" | bc -l) )); then
    echo "/gps/ene/type Gauss"
    echo "/gps/ene/mono ${energy} GeV"
    echo "/gps/ene/sigma ${sigma_E_GeV} GeV"
else
    echo "/gps/ene/type Mono"
    echo "/gps/ene/mono ${energy} GeV"
fi

echo "/run/beamOn ${nevents}"
} > ${steer_path}/${gpsmac}

# --------------------------------------------------
# STEERING FILE (PYTHON)
# --------------------------------------------------
cat > ${steer_path}/${scriptname} <<EOF
import os
from DDSim.DD4hepSimulation import DD4hepSimulation
from g4units import MeV

# Workaround for DD4hep 1.35 + Python 3.13: addParametersToRunHeader returns
# a dict that cppyy cannot convert to std::map<string,string>.
try:
    from DDSim.Helper.Meta import Meta as _M
    _M.addParametersToRunHeader = lambda self, dds: {}
except Exception:
    pass

compact_path = os.path.abspath("${geometry_folder}/SND_compact.xml")
if not os.path.isfile(compact_path):
    raise RuntimeError("Compact file not found: " + compact_path)

SIM = DD4hepSimulation()

# runType="run": simulation driven by --macroFile, not ddsim's internal gun loop.
SIM.runType   = "run"
SIM.skipNEvents = 0
SIM.compactFile  = str(compact_path)
SIM._compactFile = SIM.compactFile
SIM.outputFile   = "${sim_name}"
SIM.physicsList = "${PHYSLIST}"

# Per-chunk RNG seed: without this every chunk of a point would replay the very
# same ${nevents} events and the 50k sample would be ${nevents} events copied N times.
SIM.random.seed = ${seed}
SIM.random.enableEventSeed = True

print("COMPACT FILE  =", SIM.compactFile)
print("OUTPUT FILE   =", SIM.outputFile)
print("PARTICLE      = ${particle}")
print("Energy [GeV]  = ${energy}")
print("Beam centre   = (${pos_x}, ${pos_y}) mm")
print("sigma_x [mm]  = ${sigma_x}")
print("sigma_y [mm]  = ${sigma_y}")
print("sigma_E [frac]= ${sigma_E}")
print("chunk         = ${chunk}")
print("seed          = ${seed}")
EOF

# --------------------------------------------------
# CONDOR SHELL SCRIPT
# --------------------------------------------------
cat > ${steer_path}/${condorsh} <<EOF
#!/bin/bash
# Note: no set -e — ddsim with runType="run"+GPS crashes during Geant4 cleanup
# even when all events are written. We validate the local output file instead.

echo "Starting chunk job ${label} on \$(hostname)"

source ${repo_root}/init_key4hep.sh
export LD_LIBRARY_PATH=${repo_root}/install/lib64:${repo_root}/install/lib:\$LD_LIBRARY_PATH
export PYTHONPATH=${repo_root}/install/lib64:${repo_root}/install/lib:${repo_root}/install/python:${repo_root}:\$PYTHONPATH

LOCAL_SIM="${sim_name}"
LOCAL_TREE="${label}_ecal.root"

# ---------- 1. Geant4 simulation ----------
ddsim --enableG4GPS \\
      --macroFile    ${steer_path}/${gpsmac} \\
      --steeringFile ${steer_path}/${scriptname} \\
      &> ${log_path}/${label}.log
DDSIM_RC=\$?

if [[ ! -f "\${LOCAL_SIM}" ]] || [[ ! -s "\${LOCAL_SIM}" ]]; then
    echo "ERROR: ddsim (exit \${DDSIM_RC}) produced no usable local output file."
    exit 1
fi
if [[ \${DDSIM_RC} -ne 0 ]]; then
    echo "WARNING: ddsim exit code \${DDSIM_RC} — Geant4 cleanup crash after all events written. Local output OK."
fi

# ---------- 2. Stage out the raw sim chunk (before the fragile steps) ----------
# xrdcp, not a plain cp through the worker's /eos FUSE mount: eosxd write-back
# caching on ephemeral batch slots can report a successful local write that
# never reaches the EOS namespace before the slot is torn down.
LOCAL_SIM_SIZE=\$(stat -c%s "\${LOCAL_SIM}")
xrdcp --force "\${LOCAL_SIM}" "root://eosexperiment.cern.ch/${remote_sim}"
if [[ \$? -ne 0 ]]; then
    echo "ERROR: xrdcp stage-out of the sim chunk failed."
    exit 1
fi
REMOTE_SIM_SIZE=\$(xrdfs eosexperiment.cern.ch stat "${remote_sim}" 2>/dev/null | awk '/Size:/{print \$2}')
if [[ "\${REMOTE_SIM_SIZE}" != "\${LOCAL_SIM_SIZE}" ]]; then
    echo "ERROR: sim stage-out verification failed (local=\${LOCAL_SIM_SIZE}, remote='\${REMOTE_SIM_SIZE}')."
    exit 1
fi
echo "Raw sim chunk on EOS: ${remote_sim} (\${REMOTE_SIM_SIZE} bytes)"

# ---------- 3. Digitisation (GeV2MIP + digi + flip + channel mapping) ----------
INPUT_FILE="\${PWD}/\${LOCAL_SIM}" k4run ${repo_root}/gaudi_jobs/pid2026_common/job3_digitize.py \\
      &>> ${log_path}/${label}.log
if [[ ! -s digitized.edm4hep.root ]]; then
    echo "ERROR: digitisation produced no output (raw sim is safe on EOS, rerun the pipeline stage by hand)."
    exit 2
fi

# ---------- 4. ecal TTree conversion ----------
python3 -m analysis.sim_to_ecal_tree \\
      --input  digitized.edm4hep.root \\
      --output "\${LOCAL_TREE}" \\
      --collection SiPadHitsMapped \\
      --run ${run_number} &>> ${log_path}/${label}.log
if [[ ! -s "\${LOCAL_TREE}" ]]; then
    echo "ERROR: ecal-tree conversion produced no output (raw sim is safe on EOS)."
    exit 3
fi

# ---------- 5. Stage out the ecal chunk tree ----------
LOCAL_TREE_SIZE=\$(stat -c%s "\${LOCAL_TREE}")
xrdcp --force "\${LOCAL_TREE}" "root://eosexperiment.cern.ch/${remote_tree}"
if [[ \$? -ne 0 ]]; then
    echo "ERROR: xrdcp stage-out of the ecal tree failed."
    exit 4
fi
REMOTE_TREE_SIZE=\$(xrdfs eosexperiment.cern.ch stat "${remote_tree}" 2>/dev/null | awk '/Size:/{print \$2}')
if [[ "\${REMOTE_TREE_SIZE}" != "\${LOCAL_TREE_SIZE}" ]]; then
    echo "ERROR: ecal-tree stage-out verification failed (local=\${LOCAL_TREE_SIZE}, remote='\${REMOTE_TREE_SIZE}')."
    exit 4
fi

echo "Job finished. sim=${remote_sim}  tree=${remote_tree} (\${REMOTE_TREE_SIZE} bytes)"
EOF

chmod +x ${steer_path}/${condorsh}

# --------------------------------------------------
# CONDOR SUBMIT FILE
# --------------------------------------------------
cat > ${steer_path}/${condorsub} <<EOF
executable              = ${steer_path}/${condorsh}
log                     = ${log_path}/${condorfile}.log
output                  = ${log_path}/outfile_${condorfile}.txt
error                   = ${log_path}/errors_${condorfile}.txt
should_transfer_files   = Yes
when_to_transfer_output = ON_EXIT
# Everything big is staged out to EOS by xrdcp inside the job; do NOT let
# Condor ALSO transfer it back into the AFS steer/ dir (that duplicate is what
# filled the AFS home volume). Empty = transfer no output files back.
transfer_output_files   = ""
request_cpus            = 1
# 8 GB, not 3: a 2000-event electron chunk at 99 GeV peaks near 6 GB in the
# ecal-tree conversion (sim_to_ecal_tree materialises all frames of the chunk),
# and jobs that exceed the cgroup limit go on hold mid-campaign.
request_memory          = 8 GB
request_disk            = 12 GB
+JobFlavour             = "tomorrow"
queue 1
EOF

# --------------------------------------------------
# SUBMIT
# --------------------------------------------------
if [[ "${DRY_RUN:-0}" == "1" ]]; then
    echo "  DRY_RUN=1 — not submitting ${condorsub}"
else
    cd ${steer_path}
    condor_submit ${condorsub}
    cd - > /dev/null
fi

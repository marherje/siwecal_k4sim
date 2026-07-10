#!/bin/bash

# Crear directorios si no existen
mkdir -p steer log macros

# inputs
nevents=$1
particle=$2
energy=$3
pos_x=$4
pos_y=$5
pos_z=$6
dir_x=$7
dir_y=$8
dir_z=$9

# Paths
local=$PWD
geometry_folder="${local}/../geometry"
data_path="/eos/experiment/drdcalo/siw-ecal/TB2026-06/Simulation/Generated"
steer_path="${local}/steer"
log_path="${local}/log"

# Lista de physics lists
physl=("QGSP_BERT")

for physlist in ${physl[@]}; do

    echo "Running: particle=$particle energy=$energy position=($pos_x,$pos_y,$pos_z)"

    label=${physlist}_SND_${particle}_${energy}GeV_xyz_${pos_x}_${pos_y}_${pos_z}_dir_${dir_x}_${dir_y}_${dir_z}

    scriptname=runddsim_${label}.py
    condorsh=runddsim_${label}.sh
    condorsub=runddsim_${label}.sub
    condorfile=runddsim_${label}

    # --------------------------------------------------
    # STEERING FILE (PYTHON)
    # --------------------------------------------------
cat > ${steer_path}/${scriptname} <<EOF
import os
from DDSim.DD4hepSimulation import DD4hepSimulation
from g4units import mm, GeV

gun_direction = (${dir_x}, ${dir_y}, ${dir_z}) 
gun_position = (${pos_x} * mm, ${pos_y} * mm, ${pos_z} * mm)

compact_path = os.path.abspath("${geometry_folder}/SND_compact.xml")

if not os.path.isfile(compact_path):
    raise RuntimeError("ERROR: geometry file not found: " + compact_path)

SIM = DD4hepSimulation()

SIM.runType        = "batch"
SIM.numberOfEvents = ${nevents}
SIM.skipNEvents    = 0

SIM.compactFile = str(compact_path)
SIM._compactFile = SIM.compactFile
# Write to the job's local scratch dir first, then stage out to EOS via xrdcp
# in the condor shell script below. Writing directly to a /eos-mounted path
# from a batch worker (eosxd FUSE) can report success while the file never
# lands in the EOS namespace if the write-back cache isn't flushed before the
# job slot is torn down -- same reasoning (and fix) as generic_condor_beam.sh.
SIM.outputFile     = "output_PG_${particle}_xyz_${pos_x}_${pos_y}_${pos_z}_dir_${dir_x}_${dir_y}_${dir_z}_E${energy}.edm4hep.root"

print("COMPACT FILE =", SIM.compactFile)
print("PARTICLE =", "${particle}")
print("Energy =", ${energy})
print("Position =", gun_position)
print("Direction =", gun_direction)

SIM.enableGun      = True
SIM.gun.particle   = "${particle}"
SIM.gun.energy     = ${energy} * GeV
SIM.gun.position   = gun_position
SIM.gun.direction  = gun_direction

SIM.physicsList    = "${physlist}"
# Do NOT disable userParticleHandler: needed to write CaloHitContributions with timing.
# tracker_region_zmax/rmax are defined in the compact XML.


EOF

    # --------------------------------------------------
    # SCRIPT condor
    # --------------------------------------------------
expected_output="${data_path}/output_PG_${particle}_xyz_${pos_x}_${pos_y}_${pos_z}_dir_${dir_x}_${dir_y}_${dir_z}_E${energy}.edm4hep.root"

cat > ${steer_path}/${condorsh} <<EOF
#!/bin/bash
set -e

echo "Starting job on \$(hostname)"

source ${local}/../../init_key4hep.sh
export LD_LIBRARY_PATH=${local}/../../install/lib64:${local}/../../install/lib:\$LD_LIBRARY_PATH
export PYTHONPATH=${local}/../../install/lib64:${local}/../../install/lib:${local}/../../install/python:\$PYTHONPATH

LOCAL_OUTPUT="output_PG_${particle}_xyz_${pos_x}_${pos_y}_${pos_z}_dir_${dir_x}_${dir_y}_${dir_z}_E${energy}.edm4hep.root"
REMOTE_OUTPUT="${expected_output}"

ddsim --steeringFile ${steer_path}/${scriptname} &> ${log_path}/${label}.log

if [[ ! -s "\${LOCAL_OUTPUT}" ]]; then
    echo "ERROR: ddsim produced no usable local output file."
    exit 1
fi

# Stage out via xrdcp (synchronous xrootd transfer), NOT a plain cp through the
# worker's /eos FUSE mount -- eosxd write-back caching on ephemeral batch slots
# can report a successful write that never reaches the EOS namespace. Verify the
# remote size matches before declaring success.
LOCAL_SIZE=\$(stat -c%s "\${LOCAL_OUTPUT}")
xrdcp --force "\${LOCAL_OUTPUT}" "root://eosexperiment.cern.ch/\${REMOTE_OUTPUT}"
REMOTE_SIZE=\$(xrdfs eosexperiment.cern.ch stat "\${REMOTE_OUTPUT}" 2>/dev/null | awk '/Size:/{print \$2}')
if [[ -z "\${REMOTE_SIZE}" ]] || [[ "\${REMOTE_SIZE}" != "\${LOCAL_SIZE}" ]]; then
    echo "ERROR: stage-out verification failed (local=\${LOCAL_SIZE} bytes, remote='\${REMOTE_SIZE}')."
    exit 1
fi

echo "Job finished. Output: \${REMOTE_OUTPUT} (verified \${REMOTE_SIZE} bytes on EOS)"
EOF

    chmod +x ${steer_path}/${condorsh}

    # --------------------------------------------------
    # CONDOR SUBMIT
    # --------------------------------------------------
cat > ${steer_path}/${condorsub} <<EOF
executable              = ${steer_path}/${condorsh}
log                     = ${log_path}/${condorfile}.log
output                  = ${log_path}/outfile_${condorfile}.txt
error                   = ${log_path}/errors_${condorfile}.txt
should_transfer_files   = Yes
when_to_transfer_output = ON_EXIT
# The .edm4hep.root is staged out to EOS by xrdcp inside the job; don't let
# Condor also copy it back into the AFS steer/ dir.
transfer_output_files   = ""
+JobFlavour             = "tomorrow"
queue 1
EOF

    # --------------------------------------------------
    # SUBMIT
    # --------------------------------------------------
    cd ${steer_path}
    condor_submit ${condorsub}
    cd -

done
#!/bin/bash

# Crear directorios si no existen
mkdir -p steer data log macros

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
data_path="${local}/data"
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
SIM.outputFile     = os.path.abspath("${data_path}/output_PG_${particle}_xyz_${pos_x}_${pos_y}_${pos_z}_dir_${dir_x}_${dir_y}_${dir_z}_E${energy}.edm4hep.root")

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
cat > ${steer_path}/${condorsh} <<EOF
#!/bin/bash
set -e

echo "Starting job on \$(hostname)"

source ${local}/../../init_key4hep.sh
export LD_LIBRARY_PATH=${local}/../../install/lib64:${local}/../../install/lib:\$LD_LIBRARY_PATH
export PYTHONPATH=${local}/../../install/lib64:${local}/../../install/lib:${local}/../../install/python:\$PYTHONPATH

ddsim --steeringFile ${steer_path}/${scriptname} &> ${log_path}/${label}.log

echo "Job finished"
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
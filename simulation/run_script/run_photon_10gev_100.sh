#!/usr/bin/env bash
set -eo pipefail

repo=/home/llr/ilc/shi/code/siwecal_k4sim
data_dir=/home/llr/ilc/shi/data/siwecal_k4sim

particle=gamma
energy=10
nevents=100

pos_x=1
pos_y=1
pos_z=-1000

dir_x=0
dir_y=0
dir_z=1

steer_file="${data_dir}/steer/run_photon_10gev_100.py"
log_file="${data_dir}/log/run_photon_10gev_100.log"
output_file="${data_dir}/output/output_PG_gamma_10GeV_100evt.edm4hep.root"

mkdir -p "${data_dir}/steer" "${data_dir}/output" "${data_dir}/log"

cat > "${steer_file}" <<EOF
import os

from DDSim.DD4hepSimulation import DD4hepSimulation
from g4units import GeV, mm


repo = "${repo}"
compact_path = os.path.join(repo, "simulation", "geometry", "SND_compact.xml")
output_file = "${output_file}"

if not os.path.isfile(compact_path):
    raise RuntimeError("Geometry file not found: " + compact_path)

SIM = DD4hepSimulation()

SIM.runType = "batch"
SIM.numberOfEvents = ${nevents}
SIM.skipNEvents = 0

SIM.compactFile = compact_path
SIM._compactFile = SIM.compactFile
SIM.outputFile = output_file

SIM.enableGun = True
SIM.gun.particle = "${particle}"
SIM.gun.energy = ${energy} * GeV
SIM.gun.position = (${pos_x} * mm, ${pos_y} * mm, ${pos_z} * mm)
SIM.gun.direction = (${dir_x}, ${dir_y}, ${dir_z})

SIM.physicsList = "QGSP_BERT"
EOF

source /cvmfs/sw.hsf.org/key4hep/setup.sh -r 2026-02-01

export LD_LIBRARY_PATH="${repo}/install/lib64:${repo}/install/lib:${LD_LIBRARY_PATH:-}"
export PYTHONPATH="${repo}/install/lib64:${repo}/install/lib:${repo}/install/python:${PYTHONPATH:-}"

echo "Running foreground photon simulation"
echo "  particle=${particle}"
echo "  energy=${energy} GeV"
echo "  events=${nevents}"
echo "  output=${output_file}"
echo "  log=${log_file}"

ddsim --steeringFile "${steer_file}" > "${log_file}" 2>&1

echo "Finished: ${output_file}"

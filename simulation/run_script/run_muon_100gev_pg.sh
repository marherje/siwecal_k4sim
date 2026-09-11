#!/usr/bin/env bash
set -eo pipefail

repo=/home/llr/ilc/shi/code/siwecal_k4sim
data_dir=/home/llr/ilc/shi/data/siwecal_k4sim

particle=mu-
energy=100
angle_label=0degree
nevents=1000

# Normal incidence along +z. The start point is upstream and near the ECAL
# center, so the track stays inside the detector acceptance.
pos_x=1
pos_y=1
pos_z=-1000

dir_x=0
dir_y=0
dir_z=1

steer_dir="${data_dir}/steer/muon"
log_dir="${data_dir}/log/muon"
output_dir="${data_dir}/output/muon"

steer_file="${steer_dir}/run_muon_100gev_${angle_label}_${nevents}.py"
log_file="${log_dir}/run_muon_100gev_${angle_label}_${nevents}.log"
output_file="${output_dir}/mu-_100GeV_0degree.edm4hep.root"

mkdir -p "${steer_dir}" "${output_dir}" "${log_dir}"

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


def setup_step_limiter_apply_to_all(kernel):
    from DDG4 import PhysicsList

    physics_sequence = kernel.physicsList()
    step_limiter_physics = PhysicsList(
        kernel,
        "Geant4PhysicsList/StepLimiterApplyToAll",
    )
    step_limiter_physics.enableUI()
    physics_sequence.adopt(step_limiter_physics)

    step_limiter = step_limiter_physics.addPhysicsConstructorType(
        "G4StepLimiterPhysics"
    )
    step_limiter.SetApplyToAll(True)


SIM.physics.setupUserPhysics(setup_step_limiter_apply_to_all)
EOF

source "${repo}/init_key4hep.sh"

export LD_LIBRARY_PATH="${repo}/install/lib64:${repo}/install/lib:${LD_LIBRARY_PATH:-}"
export PYTHONPATH="${repo}/install/lib64:${repo}/install/lib:${repo}/install/python:${PYTHONPATH:-}"

echo "Running foreground muon PG simulation"
echo "  particle=${particle}"
echo "  energy=${energy} GeV"
echo "  incidence=0degree along +z"
echo "  events=${nevents}"
echo "  position=(${pos_x}, ${pos_y}, ${pos_z}) mm"
echo "  direction=(${dir_x}, ${dir_y}, ${dir_z})"
echo "  output=${output_file}"
echo "  log=${log_file}"

ddsim --steeringFile "${steer_file}" > "${log_file}" 2>&1

echo "Finished: ${output_file}"

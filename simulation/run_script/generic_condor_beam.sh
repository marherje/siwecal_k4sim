#!/bin/bash
# Beam-test simulation for SiWECAL using the Geant4 General Particle Source (GPS).
# Models a realistic test-beam: elliptical Gaussian transverse profile,
# Gaussian energy spread, and optional angular divergence.
#
# Usage:
#   ./generic_condor_beam.sh <nevents> <particle> <energy_GeV>  \
#                             <pos_x_mm> <pos_y_mm>              \
#                             <sigma_x_mm> <sigma_y_mm>          \
#                             <sigma_E_frac>                      \
#                             [theta_max_deg]
#
# Arguments:
#   nevents       Number of events per job
#   particle      Particle name (e.g. e-, mu-, pi+, proton)
#   energy_GeV    Central beam energy [GeV]
#   pos_x_mm      Beam centre X [mm]
#   pos_y_mm      Beam centre Y [mm]
#   sigma_x_mm    Gaussian sigma of the beam in X [mm]  (0 = pencil beam)
#   sigma_y_mm    Gaussian sigma of the beam in Y [mm]  (0 = pencil beam)
#   sigma_E_frac  Fractional energy spread, e.g. 0.02 = 2%  (0 = monoenergetic)
#   theta_max_deg Half-angle of flat angular divergence [deg]  (default: 0)
#
# The beam enters along +Z from z = BEAM_Z_MM.

mkdir -p steer data log

nevents=$1
particle=$2
energy=$3
pos_x=$4
pos_y=$5
sigma_x=$6
sigma_y=$7
sigma_E=$8
theta_max=${9:-0}

# Beam entry z [mm]: matches the reference simulation (10 m upstream of detector)
BEAM_Z_MM=-10000

local=$PWD
geometry_folder="${local}/../geometry"
data_path="${local}/data"
steer_path="${local}/steer"
log_path="${local}/log"

physl=("QGSP_BERT")

for physlist in ${physl[@]}; do

    echo "Beam-test run: particle=${particle}  E=${energy} GeV  centre=(${pos_x},${pos_y}) mm  sigma_x=${sigma_x} mm  sigma_y=${sigma_y} mm  sigma_E=${sigma_E}  theta_max=${theta_max} deg"

    label=${physlist}_SiWECAL_beam_${particle}_${energy}GeV_xy_${pos_x}_${pos_y}_sigx${sigma_x}_sigy${sigma_y}_sigE${sigma_E}

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
from g4units import mm, GeV, deg

compact_path = os.path.abspath("${geometry_folder}/SND_compact.xml")
if not os.path.isfile(compact_path):
    raise RuntimeError("Compact file not found: " + compact_path)

SIM = DD4hepSimulation()
SIM.runType        = "batch"
SIM.numberOfEvents = ${nevents}
SIM.skipNEvents    = 0
SIM.compactFile    = str(compact_path)
SIM._compactFile   = SIM.compactFile
SIM.outputFile     = os.path.abspath(
    "${data_path}/output_beam_${particle}_${energy}GeV_xy_${pos_x}_${pos_y}_sigx${sigma_x}_sigy${sigma_y}_sigE${sigma_E}.edm4hep.root"
)

print("COMPACT FILE  =", SIM.compactFile)
print("PARTICLE      =", "${particle}")
print("Energy [GeV]  =", ${energy})
print("Beam centre   = (${pos_x}, ${pos_y}) mm")
print("sigma_x [mm]  =", ${sigma_x})
print("sigma_y [mm]  =", ${sigma_y})
print("sigma_E [frac]=", ${sigma_E})
print("theta_max[deg]=", ${theta_max})
print("Beam z [mm]   =", ${BEAM_Z_MM})

# ── Beam simulation via Geant4 General Particle Source (GPS) ──────────────────
# GPS gives full control over the spatial, energy, and angular distributions.
# Disable the simple particle gun so only GPS fires.
SIM.enableGun = False

_sigma_E_GeV = ${energy} * ${sigma_E}   # energy sigma in absolute GeV
_theta_max   = float(${theta_max})

_gps_cmds = [
    "/gps/particle ${particle}",

    # ── Energy distribution ────────────────────────────────────────────────────
]
if ${sigma_E} > 0:
    _gps_cmds += [
        "/gps/ene/type Gauss",
        f"/gps/ene/mono ${energy} GeV",
        f"/gps/ene/sigma {_sigma_E_GeV:.6f} GeV",
    ]
else:
    _gps_cmds += [
        "/gps/ene/type Mono",
        f"/gps/ene/mono ${energy} GeV",
    ]

# ── Spatial distribution: elliptical Gaussian beam profile ────────────────────
# sigma_x and sigma_y are independent to model an asymmetric beam spot.
_gps_cmds += [
    "/gps/pos/type Beam",
    "/gps/pos/shape Circle",
    f"/gps/pos/centre ${pos_x} ${pos_y} ${BEAM_Z_MM} mm",
    f"/gps/pos/sigma_x ${sigma_x} mm",
    f"/gps/pos/sigma_y ${sigma_y} mm",
]

# ── Angular distribution: beam along +Z with optional divergence ───────────────
_gps_cmds += ["/gps/direction 0 0 1"]
if _theta_max > 0:
    # Flat (cone) distribution within theta_max around +Z.
    # For Gaussian divergence replace "beam1d" with "beam2d" and set sigma_r.
    _gps_cmds += [
        "/gps/ang/type beam1d",
        f"/gps/ang/sigma_r {_theta_max} deg",
    ]
else:
    _gps_cmds += ["/gps/ang/type beam1d"]

SIM.action.commandsConfigure = _gps_cmds

SIM.physicsList = "${physlist}"

EOF

    # --------------------------------------------------
    # CONDOR SHELL SCRIPT
    # --------------------------------------------------
cat > ${steer_path}/${condorsh} <<EOF
#!/bin/bash
set -e

echo "Starting beam-test job on \$(hostname)"

source ${local}/../../init_key4hep.sh
export LD_LIBRARY_PATH=${local}/../../install/lib64:${local}/../../install/lib:\$LD_LIBRARY_PATH
export PYTHONPATH=${local}/../../install/lib64:${local}/../../install/lib:${local}/../../install/python:\$PYTHONPATH

ddsim --steeringFile ${steer_path}/${scriptname} &> ${log_path}/${label}.log

echo "Job finished"
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

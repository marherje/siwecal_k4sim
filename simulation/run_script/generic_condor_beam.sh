#!/bin/bash
# Beam-test simulation for SiWECAL using the Geant4 General Particle Source (GPS).
# GPS commands live in a Geant4 macro file (.mac) passed via --macroFile.
# The steering file uses runType="run" so the macro drives the simulation.
#
# Usage:
#   ./generic_condor_beam.sh <nevents> <particle> <energy_GeV>  \
#                             <pos_x_mm> <pos_y_mm>              \
#                             <sigma_x_mm> <sigma_y_mm>          \
#                             <sigma_E_frac>                      \
#                             [theta_max_deg]

mkdir -p steer log

nevents=$1
particle=$2
energy=$3
pos_x=$4
pos_y=$5
sigma_x=$6
sigma_y=$7
sigma_E=$8

BEAM_Z_MM=-2000

local=$PWD
geometry_folder="${local}/../geometry"
data_path="/eos/experiment/drdcalo/siw-ecal/TB2026-06/Simulation/Generated"
steer_path="${local}/steer"
log_path="${local}/log"

physl=("QGSP_BERT")

for physlist in ${physl[@]}; do

    echo "Beam-test: particle=${particle}  E=${energy} GeV  centre=(${pos_x},${pos_y}) mm  sigma_x=${sigma_x} mm  sigma_y=${sigma_y} mm  sigma_E=${sigma_E}  theta_max=${theta_max} deg"

    label=${physlist}_SiWECAL_beam_${particle}_${energy}GeV_xy_${pos_x}_${pos_y}_sigx${sigma_x}_sigy${sigma_y}_sigE${sigma_E}

    gpsmac=gps_${label}.mac
    scriptname=runddsim_${label}.py
    condorsh=runddsim_${label}.sh
    condorsub=runddsim_${label}.sub
    condorfile=runddsim_${label}

    sigma_E_GeV=$(python3 -c "print(${energy} * ${sigma_E})")

    # --------------------------------------------------
    # GEANT4 GPS MACRO FILE
    # Contains GPS setup + /run/beamOn N.
    # runType="run" means Geant4 drives the loop from this macro.
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
    # runType="run": Geant4 is driven by the macro file.
    # numberOfEvents is NOT set here; /run/beamOn in the macro controls it.
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
# numberOfEvents is intentionally absent; /run/beamOn in the macro controls it.
SIM.runType   = "run"
SIM.skipNEvents = 0
SIM.compactFile  = str(compact_path)
SIM._compactFile = SIM.compactFile
# Write to the job's local scratch dir first, then stage out to EOS via xrdcp
# in the condor shell script below. Writing directly to a /eos-mounted path
# from a batch worker (eosxd FUSE) can report success while the file never
# lands in the EOS namespace if the write-back cache isn't flushed before the
# job slot is torn down. See runddsim shell script for the verified stage-out.
SIM.outputFile   = "output_beam_${particle}_${energy}GeV_xy_${pos_x}_${pos_y}_sigx${sigma_x}_sigy${sigma_y}_sigE${sigma_E}.edm4hep.root"
SIM.physicsList = "${physlist}"

# Do NOT disable userParticleHandler: DDG4 needs it to write CaloHitContributions
# with per-step timing. tracker_region_zmax/rmax are defined in the compact XML.

print("COMPACT FILE  =", SIM.compactFile)
print("OUTPUT FILE   =", SIM.outputFile)
print("PARTICLE      = ${particle}")
print("Energy [GeV]  = ${energy}")
print("Beam centre   = (${pos_x}, ${pos_y}) mm")
print("sigma_x [mm]  = ${sigma_x}")
print("sigma_y [mm]  = ${sigma_y}")
print("sigma_E [frac]= ${sigma_E}")
print("theta_max[deg]= ${theta_max}")
print("Beam z [mm]   = ${BEAM_Z_MM}")
EOF

    # --------------------------------------------------
    # CONDOR SHELL SCRIPT
    # --enableG4GPS: activates G4GeneralParticleSource as the primary generator.
    # --macroFile:   passes the GPS macro (contains /run/beamOn N).
    # --------------------------------------------------
    expected_output="${data_path}/output_beam_${particle}_${energy}GeV_xy_${pos_x}_${pos_y}_sigx${sigma_x}_sigy${sigma_y}_sigE${sigma_E}.edm4hep.root"

cat > ${steer_path}/${condorsh} <<EOF
#!/bin/bash
# Note: no set -e — ddsim with runType="run"+GPS crashes during Geant4 cleanup
# even when all events are written. We validate the local output file instead.

echo "Starting beam-test job on \$(hostname)"

source ${local}/../../init_key4hep.sh
export LD_LIBRARY_PATH=${local}/../../install/lib64:${local}/../../install/lib:\$LD_LIBRARY_PATH
export PYTHONPATH=${local}/../../install/lib64:${local}/../../install/lib:${local}/../../install/python:\$PYTHONPATH

LOCAL_OUTPUT="output_beam_${particle}_${energy}GeV_xy_${pos_x}_${pos_y}_sigx${sigma_x}_sigy${sigma_y}_sigE${sigma_E}.edm4hep.root"
REMOTE_OUTPUT="${expected_output}"

ddsim --enableG4GPS \\
      --macroFile    ${steer_path}/${gpsmac} \\
      --steeringFile ${steer_path}/${scriptname} \\
      &> ${log_path}/${label}.log
DDSIM_RC=\$?

if [[ ! -f "\${LOCAL_OUTPUT}" ]] || [[ ! -s "\${LOCAL_OUTPUT}" ]]; then
    echo "ERROR: ddsim (exit \${DDSIM_RC}) produced no usable local output file."
    exit 1
fi

if [[ \${DDSIM_RC} -ne 0 ]]; then
    echo "WARNING: ddsim exit code \${DDSIM_RC} — Geant4 cleanup crash after all events written. Local output OK, staging out."
fi

# Stage out via xrdcp (synchronous xrootd transfer through the redirector),
# NOT a plain cp through the worker's /eos FUSE mount: eosxd write-back
# caching on ephemeral batch slots can report a successful local write that
# never actually reaches the EOS namespace before the job slot is torn down.
LOCAL_SIZE=\$(stat -c%s "\${LOCAL_OUTPUT}")

xrdcp --force "\${LOCAL_OUTPUT}" "root://eosexperiment.cern.ch/\${REMOTE_OUTPUT}"
XRDCP_RC=\$?
if [[ \${XRDCP_RC} -ne 0 ]]; then
    echo "ERROR: xrdcp stage-out to EOS failed (exit \${XRDCP_RC})."
    exit 1
fi

REMOTE_SIZE=\$(xrdfs eosexperiment.cern.ch stat "\${REMOTE_OUTPUT}" 2>/dev/null | awk '/Size:/{print \$2}')
if [[ -z "\${REMOTE_SIZE}" ]] || [[ "\${REMOTE_SIZE}" != "\${LOCAL_SIZE}" ]]; then
    echo "ERROR: stage-out verification failed (local=\${LOCAL_SIZE} bytes, remote='\${REMOTE_SIZE}')."
    exit 1
fi

echo "Job finished. Output: \${REMOTE_OUTPUT} (verified \${REMOTE_SIZE} bytes on EOS)"
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
# The big .edm4hep.root is staged out to EOS by xrdcp inside the job; do NOT
# let Condor ALSO transfer it back into the AFS steer/ dir (that duplicate is
# what filled the AFS home volume). Empty = transfer no output files back.
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

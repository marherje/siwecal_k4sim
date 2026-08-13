#!/bin/bash
# Re-run only the post-Geant4 part of a chunk (digitisation + ACTS tracking +
# ecal-tree conversion) from a raw sim chunk that is already on EOS.
#
# Why this exists: the chunk jobs stage the raw sim to EOS *before* the
# memory-hungry conversion step, so a job that dies on the cgroup memory limit
# still leaves a perfectly good simulation behind. Releasing such a job makes
# condor redo ~50 min of Geant4 for nothing; this script skips straight to the
# cheap steps (~10 min, 8 GB).
#
# Usage:
#   ./reprocess_chunk.sh <sim_chunk_on_eos>            # submit one
#   ./reprocess_chunk.sh --missing "beam_e-_99GeV*"    # every sim chunk whose
#                                                       ecal tree is absent
#   DRY_RUN=1 ./reprocess_chunk.sh ...                 # generate, don't submit

set -uo pipefail
mkdir -p steer log

local=$PWD
repo_root="$(cd "${local}/../.." && pwd)"
eos_base="/eos/experiment/drdcalo/siw-ecal/TB2026-06/Simulation"
sim_dir="${eos_base}/Generated/chunks"
tree_dir="${eos_base}/Processed/chunks"
steer_path="${local}/steer"
log_path="${local}/log"

submit_one() {
    local simfile="$1"
    local base label condorsh condorsub
    base=$(basename "${simfile}" .edm4hep.root)   # output_<label>
    label=${base#output_}
    condorsh=reproc_${label}.sh
    condorsub=reproc_${label}.sub

    echo "Reprocessing ${label}"

cat > ${steer_path}/${condorsh} <<EOF
#!/bin/bash
echo "Reprocessing chunk ${label} on \$(hostname)"

source ${repo_root}/init_key4hep.sh
export LD_LIBRARY_PATH=${repo_root}/install/lib64:${repo_root}/install/lib:\$LD_LIBRARY_PATH
export PYTHONPATH=${repo_root}/install/lib64:${repo_root}/install/lib:${repo_root}/install/python:${repo_root}:\$PYTHONPATH

LOCAL_SIM="${base}.edm4hep.root"
LOCAL_TREE="${label}_ecal.root"
LOCAL_DIGITIZED="digitized.edm4hep.root"

xrdcp --force "root://eosexperiment.cern.ch/${sim_dir}/${base}.edm4hep.root" "\${LOCAL_SIM}"
if [[ ! -s "\${LOCAL_SIM}" ]]; then
    echo "ERROR: could not fetch the raw sim chunk from EOS."
    exit 1
fi

INPUT_FILE="\${PWD}/\${LOCAL_SIM}" k4run ${repo_root}/gaudi_jobs/pid2026_common/job3_digitize.py \\
      &> ${log_path}/reproc_${label}.log
if [[ ! -s digitized.edm4hep.root ]]; then
    echo "ERROR: digitisation produced no output."
    exit 2
fi

# Tracks on SiPadHitsDigi, BEFORE the flip: DetectorFlipper rewrites the hit z
# into the test-beam frame, which no longer matches the ACTS surfaces. Written
# to a temp file and swapped back onto digitized.edm4hep.root itself (keep *
# carries the digitised collections forward) so ACTSTracks/EMShowers/
# SiPadMeasurements end up in the ONE edm4hep file that gets staged, instead of
# a separate tracks.edm4hep.root product.
INPUT_FILE="digitized.edm4hep.root" INPUT_COLLECTION="SiPadHitsDigi" \\
      OUTPUT_FILE="\${LOCAL_DIGITIZED}.tracks_tmp" SEED_MOMENTUM=${BEAM_ENERGY} \\
      k4run ${repo_root}/gaudi_jobs/pid2026_common/job4_tracking.py \\
      &>> ${log_path}/reproc_${label}.log
if [[ ! -s "\${LOCAL_DIGITIZED}.tracks_tmp" ]]; then
    echo "ERROR: tracking produced no output."
    exit 5
fi
mv "\${LOCAL_DIGITIZED}.tracks_tmp" "\${LOCAL_DIGITIZED}"

# The run number encodes energy and chunk: E*1000 + chunk (see
# generic_condor_beam_chunk.sh), so hadd'ed trees keep chunk identity.
python3 -m analysis.sim_to_ecal_tree \\
      --input  digitized.edm4hep.root \\
      --output "\${LOCAL_TREE}" \\
      --collection SiPadHitsMapped \\
      --run ${RUN_NUMBER} &>> ${log_path}/reproc_${label}.log
if [[ ! -s "\${LOCAL_TREE}" ]]; then
    echo "ERROR: ecal-tree conversion produced no output."
    exit 3
fi

LOCAL_TREE_SIZE=\$(stat -c%s "\${LOCAL_TREE}")
xrdcp --force "\${LOCAL_TREE}" "root://eosexperiment.cern.ch/${tree_dir}/${label}_ecal.root"
REMOTE_TREE_SIZE=\$(xrdfs eosexperiment.cern.ch stat "${tree_dir}/${label}_ecal.root" 2>/dev/null | awk '/Size:/{print \$2}')
if [[ "\${REMOTE_TREE_SIZE}" != "\${LOCAL_TREE_SIZE}" ]]; then
    echo "ERROR: stage-out verification failed (local=\${LOCAL_TREE_SIZE}, remote='\${REMOTE_TREE_SIZE}')."
    exit 4
fi

LOCAL_DIGITIZED_SIZE=\$(stat -c%s "\${LOCAL_DIGITIZED}")
xrdcp --force "\${LOCAL_DIGITIZED}" "root://eosexperiment.cern.ch/${tree_dir}/${label}_digitized.edm4hep.root"
REMOTE_DIGITIZED_SIZE=\$(xrdfs eosexperiment.cern.ch stat "${tree_dir}/${label}_digitized.edm4hep.root" 2>/dev/null | awk '/Size:/{print \$2}')
if [[ "\${REMOTE_DIGITIZED_SIZE}" != "\${LOCAL_DIGITIZED_SIZE}" ]]; then
    echo "ERROR: digitized+tracks stage-out verification failed (local=\${LOCAL_DIGITIZED_SIZE}, remote='\${REMOTE_DIGITIZED_SIZE}')."
    exit 4
fi

echo "Done. tree=${tree_dir}/${label}_ecal.root (\${REMOTE_TREE_SIZE} bytes)  digitized+tracks=${tree_dir}/${label}_digitized.edm4hep.root (\${REMOTE_DIGITIZED_SIZE} bytes)"
EOF

    chmod +x ${steer_path}/${condorsh}

cat > ${steer_path}/${condorsub} <<EOF
executable              = ${steer_path}/${condorsh}
log                     = ${log_path}/reproc_${label}.condor.log
output                  = ${log_path}/outfile_reproc_${label}.txt
error                   = ${log_path}/errors_reproc_${label}.txt
should_transfer_files   = Yes
when_to_transfer_output = ON_EXIT
transfer_output_files   = ""
request_cpus            = 1
request_memory          = 8 GB
request_disk            = 20 GB
+JobFlavour             = "workday"
queue 1
EOF

    if [[ "${DRY_RUN:-0}" == "1" ]]; then
        echo "  DRY_RUN=1 — not submitting ${condorsub}"
    else
        cd ${steer_path} && condor_submit ${condorsub}; cd - > /dev/null
    fi
}

# Derive the run number (E*1000 + chunk) from the chunk label.
set_run_number() {
    local label="$1"
    local energy chunk
    energy=$(sed -E 's/.*_([0-9]+)GeV_.*/\1/' <<< "${label}")
    chunk=$(sed -E 's/.*_c([0-9]+)$/\1/' <<< "${label}")
    RUN_NUMBER=$(( energy * 1000 + 10#${chunk} ))
    # Also the seed momentum for the tracking step.
    BEAM_ENERGY=${energy}
}

if [[ "${1:-}" == "--missing" ]]; then
    pattern="${2:-*}"
    n=0
    for f in ${sim_dir}/output_${pattern}.edm4hep.root; do
        [[ -e "${f}" ]] || continue
        label=$(basename "${f}" .edm4hep.root); label=${label#output_}
        [[ -s "${tree_dir}/${label}_ecal.root" ]] && continue   # already converted
        set_run_number "${label}"
        submit_one "${f}"
        n=$(( n + 1 ))
    done
    echo "=== ${n} chunk(s) queued for reprocessing ==="
else
    [[ -n "${1:-}" ]] || { echo "Usage: $0 <sim_chunk_on_eos> | --missing <pattern>"; exit 1; }
    label=$(basename "$1" .edm4hep.root); label=${label#output_}
    set_run_number "${label}"
    submit_one "$1"
fi

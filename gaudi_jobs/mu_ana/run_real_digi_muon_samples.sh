#!/usr/bin/env bash
set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
repo_dir="$(cd "${script_dir}/../.." && pwd)"
sample_dir="/home/llr/ilc/shi/data/siwecal_k4sim/output/muon"
output_dir="${OUTPUT_DIR:-${sample_dir}/digitized}"

input_collection="${INPUT_COLLECTION:-SiPadHits}"
output_collection="${OUTPUT_COLLECTION:-SiPadHitsDigi}"
digitized_energy_collection="${DIGITIZED_ENERGY_COLLECTION:-SiPadHitsDigiDigitizedEnergy}"
digitized_time_collection="${DIGITIZED_TIME_COLLECTION:-SiPadHitsDigiDigitizedTime}"
digitized_energy_scale="${DIGITIZED_ENERGY_SCALE:-}"

mkdir -p "${output_dir}"

set +u
source /cvmfs/sw.hsf.org/key4hep/setup.sh -r 2026-02-01
if [[ -f "${repo_dir}/build/snd_simenv.sh" ]]; then
  source "${repo_dir}/build/snd_simenv.sh"
elif [[ -d "${repo_dir}/install" ]]; then
  export LD_LIBRARY_PATH="${repo_dir}/install/lib64:${repo_dir}/install/lib:${LD_LIBRARY_PATH}"
  export PYTHONPATH="${repo_dir}/install/lib64:${repo_dir}/install/lib:${repo_dir}/install/python:${PYTHONPATH}"
else
  echo "Cannot find ${repo_dir}/build/snd_simenv.sh or ${repo_dir}/install." >&2
  echo "Build the Gaudi plugin first: cd ${repo_dir} && bash build.sh" >&2
  exit 1
fi
set -u

labels=("0degree" "45degree")
inputs=(
  "${sample_dir}/mu-_100GeV_0degree.edm4hep.root"
  "${sample_dir}/mu-_100GeV_45deg.edm4hep.root"
)
outputs=(
  "${output_dir}/mu-_100GeV_0degree_real_digitized.edm4hep.root"
  "${output_dir}/mu-_100GeV_45deg_real_digitized.edm4hep.root"
)

for index in "${!labels[@]}"; do
  label="${labels[$index]}"
  input_file="${inputs[$index]}"
  output_file="${outputs[$index]}"

  if [[ ! -f "${input_file}" ]]; then
    echo "Missing input file: ${input_file}" >&2
    exit 1
  fi

  echo "Digitizing ${label} muon sample"
  echo "  input=${input_file}"
  echo "  output=${output_file}"
  if [[ -n "${digitized_energy_scale}" ]]; then
    echo "  digitized_energy_scale=${digitized_energy_scale}"
  else
    echo "  digitized_energy_scale=real_digi_muon.py default"
  fi

  run_env=(
    "INPUT_FILE=${input_file}"
    "OUTPUT_FILE=${output_file}"
    "INPUT_COLLECTION=${input_collection}"
    "OUTPUT_COLLECTION=${output_collection}"
    "DIGITIZED_ENERGY_COLLECTION=${digitized_energy_collection}"
    "DIGITIZED_TIME_COLLECTION=${digitized_time_collection}"
  )
  if [[ -n "${digitized_energy_scale}" ]]; then
    run_env+=("DIGITIZED_ENERGY_SCALE=${digitized_energy_scale}")
  fi

  env "${run_env[@]}" k4run "${script_dir}/real_digi_muon.py"
done

echo "Digitized outputs written to ${output_dir}"

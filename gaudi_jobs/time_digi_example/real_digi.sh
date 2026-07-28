#!/usr/bin/env bash
set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
repo_dir="$(cd "${script_dir}/../.." && pwd)"

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

k4run "${script_dir}/real_digi.py"

#!/usr/bin/env bash
set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
repo_dir="$(cd "${script_dir}/../.." && pwd)"

set +u
source /cvmfs/sw.hsf.org/key4hep/setup.sh -r 2026-02-01
source "${repo_dir}/build/snd_simenv.sh"
set -u

k4run "${script_dir}/real_digi.py"

#!/usr/bin/env bash
set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

INPUT_FILE="${INPUT_FILE:-/home/llr/ilc/shi/data/siwecal_k4sim/output/output_PG_gamma_10GeV_100evt.edm4hep.root}"
DIGITIZED_FILE="${DIGITIZED_FILE:-/home/llr/ilc/shi/data/siwecal_k4sim/output/output_PG_gamma_10GeV_100evt_real_digitized.edm4hep.root}"
COLLECTION="${COLLECTION:-SiPadHits}"
EVENT="${EVENT:-0}"
MAX_HITS="${MAX_HITS:-10}"

MIP_VALUE_GEV="${MIP_VALUE_GEV:-0.0002}"
MIP_THRESHOLD="${MIP_THRESHOLD:-0.1}"
DELAY_NS="${DELAY_NS:-160}"
TAU_FAST_NS="${TAU_FAST_NS:-30}"
TAU_SLOW_NS="${TAU_SLOW_NS:-180}"
ORDER_FAST="${ORDER_FAST:-2}"
ORDER_SLOW="${ORDER_SLOW:-2}"
FAST_WINDOW_NS="${FAST_WINDOW_NS:-200}"
SLOW_WINDOW_NS="${SLOW_WINDOW_NS:-500}"
FAST_NOISE_MIP="${FAST_NOISE_MIP:-0.033333333}"
SLOW_NOISE_MIP="${SLOW_NOISE_MIP:-0.083333333}"
PEAK_SEARCH_BINS="${PEAK_SEARCH_BINS:-64}"
REFINE_ITERATIONS="${REFINE_ITERATIONS:-48}"
TRIGGER_SEARCH_BINS="${TRIGGER_SEARCH_BINS:-64}"
RANDOM_SEED="${RANDOM_SEED:-5489}"

set +u
source /cvmfs/sw.hsf.org/key4hep/setup.sh -r 2026-02-01
set -u

build_dir="${script_dir}/build_key4hep"

cmake -S "${script_dir}" -B "${build_dir}" \
  -DCMAKE_PREFIX_PATH="${CMAKE_PREFIX_PATH}" \
  -DCMAKE_CXX_COMPILER="$(command -v g++)"
cmake --build "${build_dir}" -j4 --target compare_sipad_digitized

"${build_dir}/compare_sipad_digitized" \
  "${INPUT_FILE}" \
  "${DIGITIZED_FILE}" \
  "${COLLECTION}" \
  "${EVENT}" \
  "${MAX_HITS}" \
  --mip-value "${MIP_VALUE_GEV}" \
  --threshold "${MIP_THRESHOLD}" \
  --delay "${DELAY_NS}" \
  --tau-fast "${TAU_FAST_NS}" \
  --tau-slow "${TAU_SLOW_NS}" \
  --order-fast "${ORDER_FAST}" \
  --order-slow "${ORDER_SLOW}" \
  --fast-window "${FAST_WINDOW_NS}" \
  --slow-window "${SLOW_WINDOW_NS}" \
  --fast-noise "${FAST_NOISE_MIP}" \
  --slow-noise "${SLOW_NOISE_MIP}" \
  --peak-bins "${PEAK_SEARCH_BINS}" \
  --refine-iterations "${REFINE_ITERATIONS}" \
  --trigger-bins "${TRIGGER_SEARCH_BINS}" \
  --seed "${RANDOM_SEED}" \
  "$@"

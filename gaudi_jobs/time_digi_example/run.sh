#!/usr/bin/env bash
set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
repo_dir="$(cd "${script_dir}/../.." && pwd)"

# ---- Input and output -------------------------------------------------------
# Edit these defaults, or override the first four fields from the command line:
#   ./run.sh [input.root] [collection] [event] [max_hits] [extra options]
RUN_MODE="${RUN_MODE:-double}"  # single or double
INPUT_FILE="${INPUT_FILE:-/home/llr/ilc/shi/data/siwecal_k4sim/output/output_PG_mu_smoke_1evt.edm4hep.root}"
COLLECTION="${COLLECTION:-SiPadHits}"
EVENT="${EVENT:-0}"
MAX_HITS="${MAX_HITS:-1}"  # only used by RUN_MODE=single
OUTPUT_DIR="${OUTPUT_DIR:-${script_dir}/figures}"
CLEAN_FIGURES="${CLEAN_FIGURES:-1}"
PARTICLE="${PARTICLE:-photon}"
ENERGY_GEV="${ENERGY_GEV:-10}"
SINGLE_OUTPUT_STEM="${SINGLE_OUTPUT_STEM:-${PARTICLE}_${ENERGY_GEV}GeV_hit}"
DOUBLE_OUTPUT_STEM="${DOUBLE_OUTPUT_STEM:-${PARTICLE}_${ENERGY_GEV}GeV_double_hit}"

# ---- Shaping configuration -------------------------------------------------
# These values are passed to plot_sipad_shaping and mirror CellShapingConfig.
MIP_VALUE_GEV=0.0002
MIP_THRESHOLD=0.5
DELAY_NS=200
TAU_FAST_NS=50
TAU_SLOW_NS=180
ORDER_FAST=2
ORDER_SLOW=2
FAST_WINDOW_NS=200
SLOW_WINDOW_NS=500
FAST_NOISE_MIP=0.033333333
SLOW_NOISE_MIP=0.083333333
PEAK_SEARCH_BINS=64
REFINE_ITERATIONS=48
TRIGGER_SEARCH_BINS=64
RANDOM_SEED=5489

if (($# >= 1)); then INPUT_FILE="$1"; fi
if (($# >= 2)); then COLLECTION="$2"; fi
if (($# >= 3)); then EVENT="$3"; fi
if (($# >= 4)); then MAX_HITS="$4"; fi
if (($# >= 5)); then
  EXTRA_ARGS=("${@:5}")
else
  EXTRA_ARGS=()
fi

set +u
source /cvmfs/sw.hsf.org/key4hep/setup.sh -r 2026-02-01
set -u

build_dir="${script_dir}/build_key4hep"

case "${RUN_MODE}" in
  single)
    executable="${build_dir}/plot_sipad_shaping"
    output_stem="${SINGLE_OUTPUT_STEM}"
    ;;
  double)
    executable="${build_dir}/doublehit"
    output_stem="${DOUBLE_OUTPUT_STEM}"
    ;;
  *)
    echo "Unknown RUN_MODE='${RUN_MODE}'. Use 'single' or 'double'." >&2
    exit 1
    ;;
esac

mkdir -p "${OUTPUT_DIR}"
if [[ "${CLEAN_FIGURES}" == "1" ]]; then
  find "${OUTPUT_DIR}" -maxdepth 1 \
    \( -name 'hit_*.pdf' \
    -o -name 'merged_hit_*.pdf' \
    -o -name 'event_*_contribution_time_vs_z.pdf' \
    -o -name "${output_stem}.pdf" \
    -o -name "${output_stem}_*.pdf" \
    -o -name 'shaping_summary.csv' \) -delete
fi

cmake -S "${script_dir}" -B "${build_dir}" \
  -DCMAKE_PREFIX_PATH="${CMAKE_PREFIX_PATH}" \
  -DCMAKE_CXX_COMPILER="$(command -v g++)"
cmake --build "${build_dir}" -j4

common_args=(
  --mip-value "${MIP_VALUE_GEV}"
  --threshold "${MIP_THRESHOLD}"
  --delay "${DELAY_NS}"
  --tau-fast "${TAU_FAST_NS}"
  --tau-slow "${TAU_SLOW_NS}"
  --order-fast "${ORDER_FAST}"
  --order-slow "${ORDER_SLOW}"
  --fast-window "${FAST_WINDOW_NS}"
  --slow-window "${SLOW_WINDOW_NS}"
  --fast-noise "${FAST_NOISE_MIP}"
  --slow-noise "${SLOW_NOISE_MIP}"
  --peak-bins "${PEAK_SEARCH_BINS}"
  --refine-iterations "${REFINE_ITERATIONS}"
  --trigger-bins "${TRIGGER_SEARCH_BINS}"
  --seed "${RANDOM_SEED}"
  --output-stem "${output_stem}"
)

if [[ "${RUN_MODE}" == "single" ]]; then
  "${executable}" \
    "${INPUT_FILE}" \
    "${OUTPUT_DIR}" \
    "${COLLECTION}" \
    "${EVENT}" \
    "${MAX_HITS}" \
    "${common_args[@]}" \
    "${EXTRA_ARGS[@]}"
else
  "${executable}" \
    "${INPUT_FILE}" \
    "${OUTPUT_DIR}" \
    "${COLLECTION}" \
    "${EVENT}" \
    "${common_args[@]}" \
    "${EXTRA_ARGS[@]}"
fi

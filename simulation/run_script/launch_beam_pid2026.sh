#!/bin/bash
# Launch the 2026 PID campaign: e-, mu- and pi- beams at 74 and 99 GeV,
# beam centre (1.139, 1.164) mm, 50 000 events per (particle, energy, position).
#
# Each point is split into NCHUNKS condor jobs of CHUNK_EVENTS events, each with
# its own RNG seed (a 50k single job would be ~26 h of wall time for electrons).
# Every job runs ddsim + digitisation + ecal-tree conversion and stages out
#   Generated/chunks/output_<point>_cNNN.edm4hep.root
#   Processed/chunks/<point>_cNNN_ecal.root
#
# Usage:
#   source init_siwecal_soft.sh
#   cd simulation/run_script
#   bash launch_beam_pid2026.sh                 # submit everything
#   DRY_RUN=1 bash launch_beam_pid2026.sh       # generate steer files only
#   PARTICLES="pi-" ENERGIES="99" bash launch_beam_pid2026.sh   # subset / resubmit
#   CHUNKS="3 17" bash launch_beam_pid2026.sh   # re-run individual failed chunks

set -uo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")"

TOTAL_EVENTS=${TOTAL_EVENTS:-50000}
CHUNK_EVENTS=${CHUNK_EVENTS:-2000}
NCHUNKS=$(( TOTAL_EVENTS / CHUNK_EVENTS ))

PARTICLES=${PARTICLES:-"e- mu- pi-"}
ENERGIES=${ENERGIES:-"74 99"}
CHUNKS=${CHUNKS:-$(seq 0 $(( NCHUNKS - 1 )))}

pos_x=1.139   # mm — beam centre x
pos_y=1.164   # mm — beam centre y

sigma_E=0.02  # fractional energy spread (2 %), same for all species

# Per-species beam optics (mm), as used in launch_beam.sh
sigma_x_e=13.75;  sigma_y_e=8.25
sigma_x_mu=38.5;  sigma_y_mu=46.75
sigma_x_pi=24.75; sigma_y_pi=13.75

# Deterministic, collision-free seeds: species offset + energy + chunk.
seed_offset_for() {
    case "$1" in
        "e-")  echo 1000000 ;;
        "mu-") echo 2000000 ;;
        "pi-") echo 3000000 ;;
        *)     echo 9000000 ;;
    esac
}

nsub=0
for particle in ${PARTICLES}; do
    case "${particle}" in
        "e-")  sigma_x=${sigma_x_e};  sigma_y=${sigma_y_e}  ;;
        "mu-") sigma_x=${sigma_x_mu}; sigma_y=${sigma_y_mu} ;;
        "pi-") sigma_x=${sigma_x_pi}; sigma_y=${sigma_y_pi} ;;
        *) echo "Unknown particle '${particle}'"; exit 1 ;;
    esac
    offset=$(seed_offset_for "${particle}")

    for energy in ${ENERGIES}; do
        echo "=== ${particle} @ ${energy} GeV — ${NCHUNKS} chunks x ${CHUNK_EVENTS} events at (${pos_x}, ${pos_y}) mm ==="
        for chunk in ${CHUNKS}; do
            seed=$(( offset + energy * 1000 + chunk ))
            ./generic_condor_beam_chunk.sh \
                ${CHUNK_EVENTS} ${particle} ${energy} \
                ${pos_x} ${pos_y} \
                ${sigma_x} ${sigma_y} \
                ${sigma_E} ${chunk} ${seed} || echo "  !! submit failed for chunk ${chunk}"
            nsub=$(( nsub + 1 ))
        done
    done
done

echo
echo "=== ${nsub} chunk jobs handled (${TOTAL_EVENTS} events per point) ==="
echo "Monitor with: condor_q \$USER"

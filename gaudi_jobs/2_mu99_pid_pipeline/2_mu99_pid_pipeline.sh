#!/bin/bash
# PID 2026 pipeline — mu- beam at 99 GeV, beam centre (1.139, 1.164) mm, 50k events/energy.
#
# Inputs : Processed/chunks/*_ecal.root produced by the condor chunk jobs
#          (simulation/run_script/launch_beam_pid2026.sh)
# Outputs: Processed/<label>_ecal.{root,edm4hep.root,valtree.root}
#
# Usage:  bash gaudi_jobs/2_mu99_pid_pipeline/2_mu99_pid_pipeline.sh [--allow-partial]
set -uo pipefail
exec bash "$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)/pid2026_common/pipeline_common.sh" \
    --particle mu- --energy 99 "$@"

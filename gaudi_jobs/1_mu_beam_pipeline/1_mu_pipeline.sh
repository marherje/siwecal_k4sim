#!/bin/bash
# Full siwecal_k4sim muon beam pipeline.
#
# Steps 1-5: simulation → digitization → tracking → RNTuple (existing)
# Step 6:    digitized hits → ecal TTree → k4SiWEcalReco (shower variables)
#
# Run from the repo root:
#   source init_key4hep.sh
#   bash gaudi_jobs/1_mu_beam_pipeline/1_mu_pipeline.sh
#
# After step 6, launch the event viewer with:
#   bash event_viewer/launch_sim.sh

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
cd "${REPO_ROOT}/gaudi_jobs/1_mu_beam_pipeline"

# --- Steps 1-5: existing Gaudi jobs ---
k4run job1_shuffler.py
k4run job2_splitter.py
INPUT_FILE=timewindows.edm4hep.root k4run job3_digitize.py
k4run job4_tracking.py
k4run job5_rntuple.py

# --- Step 6: shower variables via k4SiWEcalReco ---
cd "${REPO_ROOT}"
bash analysis/run_pid_sim.sh --format both

#!/bin/bash
cd /afs/cern.ch/user/m/marquezh/public/siwecal_k4sim && source init_key4hep.sh >/dev/null 2>&1 && python3 analysis/cm_sag_collect.py /eos/experiment/drdcalo/siw-ecal/TB2026-06/Simulation/tests/sim_tuning/occupancy/sag2 --run $1 --chunks 10-17

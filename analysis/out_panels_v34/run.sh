#!/bin/bash
cd /afs/cern.ch/user/m/marquezh/public/siwecal_k4sim
source /afs/cern.ch/user/m/marquezh/public/siwecal-tb2026/setupkey4hep.sh >/dev/null 2>&1
nice python3 analysis/final_panels.py analysis/out_panels_v34
echo "EXIT=$?"

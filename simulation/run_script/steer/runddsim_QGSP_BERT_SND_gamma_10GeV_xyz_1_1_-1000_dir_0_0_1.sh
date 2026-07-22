#!/bin/bash
set -e

echo "Starting job on $(hostname)"

source /home/llr/ilc/shi/code/siwecal_k4sim/simulation/run_script/../../init_key4hep.sh
export LD_LIBRARY_PATH=/home/llr/ilc/shi/code/siwecal_k4sim/simulation/run_script/../../install/lib64:/home/llr/ilc/shi/code/siwecal_k4sim/simulation/run_script/../../install/lib:$LD_LIBRARY_PATH
export PYTHONPATH=/home/llr/ilc/shi/code/siwecal_k4sim/simulation/run_script/../../install/lib64:/home/llr/ilc/shi/code/siwecal_k4sim/simulation/run_script/../../install/lib:/home/llr/ilc/shi/code/siwecal_k4sim/simulation/run_script/../../install/python:$PYTHONPATH

ddsim --steeringFile /home/llr/ilc/shi/code/siwecal_k4sim/simulation/run_script/steer/runddsim_QGSP_BERT_SND_gamma_10GeV_xyz_1_1_-1000_dir_0_0_1.py &> /home/llr/ilc/shi/code/siwecal_k4sim/simulation/run_script/log/QGSP_BERT_SND_gamma_10GeV_xyz_1_1_-1000_dir_0_0_1.log

echo "Job finished"

#!/usr/bin/env bash

# To use only if the project is already compile and you want to play with the steering file, for example. 
# If you want to compile the project, please follow the instructions in the README.md file.
#  -r 2026-02-01
source /cvmfs/sw.hsf.org/key4hep/setup.sh -r 2026-02-01
export LD_LIBRARY_PATH=$PWD/install/lib64:$PWD/install/lib:$LD_LIBRARY_PATH
export PYTHONPATH=$PWD/install/lib64:$PWD/install/lib:$PWD/install/python:$PYTHONPATH

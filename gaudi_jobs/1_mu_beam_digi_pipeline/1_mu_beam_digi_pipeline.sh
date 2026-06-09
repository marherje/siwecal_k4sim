#!/bin/bash
set -e
k4run job1_shuffler.py
k4run job2_splitter.py
INPUT_FILE=timewindows.edm4hep.root k4run job3_digitize.py
k4run job4_tracking.py
k4run job5_rntuple.py

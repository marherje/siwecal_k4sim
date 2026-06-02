#!/bin/bash
k4run job1_shuffler.py
k4run job2_splitter.py
sleep 1
INPUT_FILE=timewindows.edm4hep.root k4run job3_digitize.py
sleep 1
k4run job4_tracking.py
sleep 1
k4run job5_rntuple.py
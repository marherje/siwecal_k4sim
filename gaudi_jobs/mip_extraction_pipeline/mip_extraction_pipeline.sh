#!/bin/bash
# MIP extraction pipeline.
# Step 1-2: shuffle and split simulation output.
# Step 3:   run MIPExtractor → produces mip_values.txt + mip_extraction.root.
# The extracted values in mip_values.txt can be pasted into job3_digitize.py
# of the reconstruction pipeline as MIPValues = [...].
set -e
k4run job1_shuffler.py
k4run job2_splitter.py
k4run job3_mipextract.py
echo "MIP extraction complete. Results in mip_values.txt and mip_extraction.root"

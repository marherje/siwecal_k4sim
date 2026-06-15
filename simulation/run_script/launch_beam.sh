#!/bin/bash
# Launch beam-test simulation jobs on HTCondor.
# Mirrors launch_PG.sh but uses the GPS beam setup (generic_condor_beam.sh).
#
# Beam parameters (sigma_x, sigma_y) match typical CERN SPS H8 beamline optics.

nevents=1000

sigma_x=23.5   # mm — horizontal beam sigma
sigma_y=29.7   # mm — vertical beam sigma
sigma_E=0.02   # fractional energy spread (2 %)
theta_max=0    # deg — no angular divergence

for particle in "mu-" "e-"
do
    for energy in 5
    do
        for pos_x in 1
        do
            for pos_y in 1
            do
                echo "Submitting: ${particle} ${energy} GeV at (${pos_x}, ${pos_y}) mm"
                ./generic_condor_beam.sh \
                    $nevents $particle $energy \
                    $pos_x $pos_y \
                    $sigma_x $sigma_y \
                    $sigma_E $theta_max
            done
        done
    done
done

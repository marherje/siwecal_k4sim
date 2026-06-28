#!/bin/bash
# Launch beam-test simulation jobs on HTCondor.
# Mirrors launch_PG.sh but uses the GPS beam setup (generic_condor_beam.sh).
#
# Beam parameters (sigma_x, sigma_y) match typical CERN SPS H8 beamline optics.

nevents=1000

sigma_x_e=13.75   # mm — horizontal beam sigma
sigma_y_e=8.25  # mm — vertical beam sigma
sigma_E_e=0.02   # fractional energy spread (2 %)

sigma_x_mu=20.5   # mm — horizontal beam sigma
sigma_y_mu=16.5   # mm — vertical beam sigma
sigma_E_mu=0.02   # fractional energy spread (2 %)

theta_max=0    # deg — no angular divergence



for particle in "e-"
do
    for energy in 54
    do
        for pos_x in -45
        do
            for pos_y in 45
            do
                echo "Submitting: ${particle} ${energy} GeV at (${pos_x}, ${pos_y}) mm"
                ./generic_condor_beam.sh \
                    $nevents $particle $energy \
                    $pos_x $pos_y \
                    $sigma_x_e $sigma_y_e \
                    $sigma_E_e $theta_max
            done
        done
    done
done


for particle in "mu-"
do
    for energy in 100
    do
        for pos_x in 1
        do
            for pos_y in 1
            do
                echo "Submitting: ${particle} ${energy} GeV at (${pos_x}, ${pos_y}) mm"
                ./generic_condor_beam.sh \
                    $nevents $particle $energy \
                    $pos_x $pos_y \
                    $sigma_x_mu $sigma_y_mu \
                    $sigma_E_mu $theta_max
            done
        done
    done
done

#!/bin/bash

nevents=100
pos_z=-1000

dir_x=0

for particle in "mu-" #"e-" #"kaon-" #"e-" "mu-" "pi-"
do
    for energy in 5 50 #4 6 8 10 20 30 40 50 60 70 80 90 100 125 150 175 200
    do
        for pos_x in 1 82.5 #0.1 0.5 1.0 2.0 5.0 10.0
        do
            for pos_y in 1 82.5 -82.5 #0.1 0.5 1.0 2.0 5.0 10.0
            do
                for dir_z in 1
                do
                    for dir_y in 0 0.05
                    do
                        echo "Running simulation for particle ${particle} with energy ${energy} GeV at position (${pos_x}, ${pos_y}, ${pos_z}) mm"
                        echo "Direction: (${dir_x}, ${dir_y}, ${dir_z})"
                        ./generic_condor.sh $nevents $particle $energy $pos_x $pos_y $pos_z $dir_x $dir_y $dir_z
                    done
                done
            done
        done
    done
done

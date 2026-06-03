#!/bin/bash
# Lanza simulaciones de beam test usando generic_condor_beam.sh (GPS de Geant4).
# El beam entra siempre a lo largo de +Z; su posición z está fija en el script.

nevents=100

# Perfil del beam — valores de la simulación de referencia
sigma_x=23.5   # sigma del beam en X [mm]
sigma_y=29.7   # sigma del beam en Y [mm]
sigma_E=0.02   # dispersión de energía (fracción), e.g. 0.02 = 2 %
theta_max=0    # divergencia angular máxima [deg]; 0 = beam perfectamente paralelo
               # Descomentar para simular divergencia: theta_max=0.1

for particle in "mu-" #"e-" #"kaon-" #"e-" "mu-" "pi-"
do
    for energy in 5 50 #4 6 8 10 20 30 40 50 60 70 80 90 100 125 150 175 200
    do
        for pos_x in 1 82.5 #0.1 0.5 1.0 2.0 5.0 10.0
        do
            for pos_y in 1 82.5 -82.5 #0.1 0.5 1.0 2.0 5.0 10.0
            do
                echo "Beam test: ${particle}  E=${energy} GeV  centre=(${pos_x}, ${pos_y}) mm  sigma_x=${sigma_x} mm  sigma_y=${sigma_y} mm  sigma_E=${sigma_E}  theta_max=${theta_max} deg"
                ./generic_condor_beam.sh $nevents $particle $energy $pos_x $pos_y $sigma_x $sigma_y $sigma_E $theta_max
            done
        done
    done
done

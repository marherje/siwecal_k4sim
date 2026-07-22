import os
from DDSim.DD4hepSimulation import DD4hepSimulation
from g4units import mm, GeV

gun_direction = (0, 0, 1) 
gun_position = (1 * mm, 1 * mm, -1000 * mm)

compact_path = os.path.abspath("/home/llr/ilc/shi/code/siwecal_k4sim/simulation/run_script/../geometry/SND_compact.xml")

if not os.path.isfile(compact_path):
    raise RuntimeError("ERROR: geometry file not found: " + compact_path)

SIM = DD4hepSimulation()

SIM.runType        = "batch"
SIM.numberOfEvents = 10
SIM.skipNEvents    = 0

SIM.compactFile = str(compact_path)
SIM._compactFile = SIM.compactFile
SIM.outputFile     = os.path.abspath("/home/llr/ilc/shi/code/siwecal_k4sim/simulation/run_script/data/output_PG_gamma_xyz_1_1_-1000_dir_0_0_1_E10.edm4hep.root")

print("COMPACT FILE =", SIM.compactFile)
print("PARTICLE =", "gamma")
print("Energy =", 10)
print("Position =", gun_position)
print("Direction =", gun_direction)

SIM.enableGun      = True
SIM.gun.particle   = "gamma"
SIM.gun.energy     = 10 * GeV
SIM.gun.position   = gun_position
SIM.gun.direction  = gun_direction

SIM.physicsList    = "QGSP_BERT"
# Do NOT disable userParticleHandler: needed to write CaloHitContributions with timing.
# tracker_region_zmax/rmax are defined in the compact XML.



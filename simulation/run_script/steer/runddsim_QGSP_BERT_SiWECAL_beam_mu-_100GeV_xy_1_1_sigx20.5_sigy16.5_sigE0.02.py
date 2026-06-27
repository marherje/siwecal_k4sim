import os
from DDSim.DD4hepSimulation import DD4hepSimulation
from g4units import MeV

# Workaround for DD4hep 1.35 + Python 3.13: addParametersToRunHeader returns
# a dict that cppyy cannot convert to std::map<string,string>.
try:
    from DDSim.Helper.Meta import Meta as _M
    _M.addParametersToRunHeader = lambda self, dds: {}
except Exception:
    pass

compact_path = os.path.abspath("/afs/cern.ch/user/m/marquezh/public/siwecal_k4sim/simulation/run_script/../geometry/SND_compact.xml")
if not os.path.isfile(compact_path):
    raise RuntimeError("Compact file not found: " + compact_path)

SIM = DD4hepSimulation()

# runType="run": simulation driven by --macroFile, not ddsim's internal gun loop.
# numberOfEvents is intentionally absent; /run/beamOn in the macro controls it.
SIM.runType   = "run"
SIM.skipNEvents = 0
SIM.compactFile  = str(compact_path)
SIM._compactFile = SIM.compactFile
SIM.outputFile   = os.path.abspath(
    "/afs/cern.ch/user/m/marquezh/public/siwecal_k4sim/simulation/run_script/data/output_beam_mu-_100GeV_xy_1_1_sigx20.5_sigy16.5_sigE0.02.edm4hep.root"
)
SIM.physicsList = "QGSP_BERT"

# Do NOT disable userParticleHandler: DDG4 needs it to write CaloHitContributions
# with per-step timing. tracker_region_zmax/rmax are defined in the compact XML.

print("COMPACT FILE  =", SIM.compactFile)
print("OUTPUT FILE   =", SIM.outputFile)
print("PARTICLE      = mu-")
print("Energy [GeV]  = 100")
print("Beam centre   = (1, 1) mm")
print("sigma_x [mm]  = 20.5")
print("sigma_y [mm]  = 16.5")
print("sigma_E [frac]= 0.02")
print("theta_max[deg]= ")
print("Beam z [mm]   = -2000")

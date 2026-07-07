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
# Write to the job's local scratch dir first, then stage out to EOS via xrdcp
# in the condor shell script below. Writing directly to a /eos-mounted path
# from a batch worker (eosxd FUSE) can report success while the file never
# lands in the EOS namespace if the write-back cache isn't flushed before the
# job slot is torn down. See runddsim shell script for the verified stage-out.
SIM.outputFile   = "output_beam_e-_74GeV_xy_-45_45_sigx13.75_sigy8.25_sigE0.02.edm4hep.root"
SIM.physicsList = "QGSP_BERT"

# Do NOT disable userParticleHandler: DDG4 needs it to write CaloHitContributions
# with per-step timing. tracker_region_zmax/rmax are defined in the compact XML.

print("COMPACT FILE  =", SIM.compactFile)
print("OUTPUT FILE   =", SIM.outputFile)
print("PARTICLE      = e-")
print("Energy [GeV]  = 74")
print("Beam centre   = (-45, 45) mm")
print("sigma_x [mm]  = 13.75")
print("sigma_y [mm]  = 8.25")
print("sigma_E [frac]= 0.02")
print("theta_max[deg]= ")
print("Beam z [mm]   = -2000")

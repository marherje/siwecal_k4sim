from k4FWCore import ApplicationMgr, IOSvc
from Configurables import ContributionInspector
from Gaudi.Configuration import DEBUG

# Input file -- change this to point to any simulation output
INPUT_FILE = "output.edm4hep.root"

iosvc = IOSvc()
iosvc.Input = INPUT_FILE
# No output file: this algorithm is purely diagnostic

inspector = ContributionInspector("ContribInspector")
inspector.SiTargetCollection = "SiPadHits"
inspector.SiPadCollection    = "SiPadHits"
inspector.MaxHitsToPrint     = 20
inspector.OutputLevel        = DEBUG  # show per-contribution DEBUG lines

ApplicationMgr(
    EvtSel = "NONE",
    EvtMax = 10,
    TopAlg = [inspector],
    ExtSvc = [iosvc],
)

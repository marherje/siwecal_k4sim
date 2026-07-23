from k4FWCore import ApplicationMgr, IOSvc
from Configurables import MIPExtractor
from Gaudi.Configuration import INFO

SIPAD_BITFIELD = "system:8,layer:8,slice:5,x:9,y:9"
SIPAD_NLAYERS  = 15

iosvc = IOSvc()
iosvc.Input = "timewindows.edm4hep.root"

extractor = MIPExtractor("MIPExtractor")
extractor.InputCollection = "SiPadHitsWindowed"  # raw hits in GeV (before GeV2MIP)
extractor.BitField        = SIPAD_BITFIELD
extractor.NLayers         = SIPAD_NLAYERS
extractor.FitMode         = 1         # 1=Gaussian, 2=Landau, 3=LanGaus (best)
extractor.MinEntries      = 100        # minimum hits per layer to attempt fit
extractor.OutputRootFile  = "mip_extraction.root"
extractor.OutputTextFile  = "mip_values.txt"
extractor.OutputLevel     = INFO

ApplicationMgr(
    EvtSel  = "NONE",
    EvtMax  = -1,
    TopAlg  = [extractor],
    ExtSvc  = [iosvc]
)

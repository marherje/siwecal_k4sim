from k4FWCore import ApplicationMgr
from Configurables import EventWindowSplitter

splitter = EventWindowSplitter("EventWindowSplitter")
splitter.InputFile             = "shuffled.edm4hep.root"
splitter.InputCollectionSiPad  = "SiPadHitsMerged"
splitter.OutputFile            = "timewindows.edm4hep.root"
splitter.OutputCollectionSiPad = "SiPadHitsWindowed"
splitter.WindowSize = 30.0   # ns

ApplicationMgr(
    EvtSel  = "NONE",
    EvtMax  = 1,
    TopAlg  = [splitter],
    ExtSvc  = []
)

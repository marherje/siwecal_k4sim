from k4FWCore import ApplicationMgr
from Configurables import EventShuffler

shuffler = EventShuffler("EventShuffler")
shuffler.InputFiles = [
    "/eos/experiment/drdcalo/siw-ecal/TB2026-06/Simulation/Generated/output_beam_e-_54GeV_xy_-45_45_sigx13.75_sigy8.25_sigE0.02.edm4hep.root",
]
shuffler.SourceIDs = [1]
shuffler.Delays    = [30.0000000001]   # ns
shuffler.CollectionsSiPad      = ["SiPadHits"]
shuffler.MaxEventsPerSource    = 1000
shuffler.OutputFile            = "shuffled.edm4hep.root"
shuffler.OutputCollectionSiPad = "SiPadHitsMerged"

ApplicationMgr(
    EvtSel  = "NONE",
    EvtMax  = 1,
    TopAlg  = [shuffler],
    ExtSvc  = []
)

from k4FWCore import ApplicationMgr
from Configurables import EventShuffler

# EventShuffler reads input files directly via podio::ROOTReader.
# IOSvc is NOT used here: do not set iosvc.Input.
# All work happens in finalize() after Gaudi processes exactly 1 dummy event.

shuffler = EventShuffler("EventShuffler")
shuffler.InputFiles = [
    "/eos/experiment/drdcalo/siw-ecal/TB2026-06/Simulation/Generated/output_mu-_E5.edm4hep.root",
]
shuffler.SourceIDs = [1]
shuffler.Delays    = [30.0000001]   # ns
shuffler.CollectionsSiPad      = ["SiPadHits"]
shuffler.MaxEventsPerSource    = 50
shuffler.OutputFile            = "shuffled.edm4hep.root"
shuffler.OutputCollectionSiPad = "SiPadHitsMerged"

ApplicationMgr(
    EvtSel  = "NONE",
    EvtMax  = 1,
    TopAlg  = [shuffler],
    ExtSvc  = []
)

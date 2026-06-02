from k4FWCore import ApplicationMgr
from Configurables import EventShuffler

# EventShuffler reads input files directly via podio::ROOTReader.
# IOSvc is NOT used here: do not set iosvc.Input.
# All work happens in finalize() after Gaudi processes exactly 1 dummy event.

shuffler = EventShuffler("EventShuffler")
shuffler.InputFiles = [
    "../../simulation/run_script/data/output_mu-_xyz_1_1_-1000_dir_0_0.05_1_E50.edm4hep.root",
    "../../simulation/run_script/data/output_mu-_xyz_82.5_82.5_-1000_dir_0_0_1_E50.edm4hep.root"
]
shuffler.SourceIDs = [1, 2]
shuffler.Delays    = [25.0001, 15.0]   # ns, adjust per source
shuffler.CollectionsSiPad = [
    "SiPadHits",
    "SiPadHits"
]
shuffler.MaxEventsPerSource       = 50
shuffler.OutputFile               = "shuffled.edm4hep.root"
shuffler.OutputCollectionSiPad  = "SiPadHitsMerged"

ApplicationMgr(
    EvtSel  = "NONE",
    EvtMax  = 1,       # Only 1 Gaudi event: execute() is a no-op, all work is in finalize()
    TopAlg  = [shuffler],
    ExtSvc  = []
)

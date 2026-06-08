from k4FWCore import ApplicationMgr
from Configurables import EventShuffler

# EventShuffler reads input files directly via podio::ROOTReader.
# IOSvc is NOT used here. All work happens in finalize() after 1 dummy event.

shuffler = EventShuffler("EventShuffler")
shuffler.InputFiles = [
    "../../simulation/run_script/data/output_beam_mu-_50GeV_xy_1_1_sigx23.5_sigy29.7_sigE0.02.edm4hep.root",
]
shuffler.SourceIDs = [1]
shuffler.Delays    = [30.0000001]   # ns
shuffler.CollectionsSiPad      = ["SiPadHits"]
shuffler.MaxEventsPerSource    = 200
shuffler.OutputFile            = "shuffled.edm4hep.root"
shuffler.OutputCollectionSiPad = "SiPadHitsMerged"

ApplicationMgr(
    EvtSel  = "NONE",
    EvtMax  = 1,
    TopAlg  = [shuffler],
    ExtSvc  = []
)

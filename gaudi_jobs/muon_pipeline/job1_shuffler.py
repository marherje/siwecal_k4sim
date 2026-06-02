from k4FWCore import ApplicationMgr
from Configurables import EventShuffler

# EventShuffler reads input files directly via podio::ROOTReader.
# IOSvc is NOT used here: do not set iosvc.Input.
# All work happens in finalize() after Gaudi processes exactly 1 dummy event.

shuffler = EventShuffler("EventShuffler")
shuffler.InputFiles = [
    # "../../output_mu-_xyz_1_-42.5_-1000_dir_0_0_1_E75.edm4hep.root",
    # "../../output_mu-_xyz_1_-42.5_-1000_dir_0_0_1_E75_1000events.edm4hep.root"
    # "../../output_e-_xyz_1_-82.5_-1000_dir_0_0_1_E75.edm4hep.root"
    # "../../output_mu-_xyz_1_-42.5_-1000_dir_0_0_1_E30_100events.edm4hep.root"
    # "../../output_mu-_xyz_1_-42.5_-1000_dir_0_0_1_E5_100events.edm4hep.root"
    # "../../output/output_mu-_xyz_1_-42.5_750_dir_0_0_1_E5_100only_mtc_events.edm4hep.root"
    # "../../output/output_mu-_xyz_1.0_-82.5_-1000.0_dir_0_0_1_E5_100_events.edm4hep.root"
    # "../../output/output_mu-_xyz_1.0_-0.0_-1000.0_dir_0_0_1_E5_100_events.edm4hep.root"
    # "../../output/output_mu-_xyz_5.0_-82.5_750.0_dir_0_0_1_E5_100_events.edm4hep.root" # working scifi hits
    # "../../output/output_mu-_xyz_5.0_-82.5_750.0_dir_0_0_1_E5_100_events.edm4hep.root" # no magfield
    # "../../output/output_mu-_xyz_5.0_-82.5_-1000.0_dir_0_0_1_E5_100_events.edm4hep.root"
    # "../../output/output_mu-_xyz_5.0_-82.5_275.0_dir_0_0_1_E5_100_events.edm4hep.root"
    #"../../output/output_mu-_xyz_5.0_-82.5_-900.0_dir_0_0_1_E5_100_events.edm4hep.root"
    "../../simulation/run_script/data/output_mu-_xyz_1_-82.5_-1000_dir_0_0_1_E5.edm4hep.root"
]
shuffler.SourceIDs = [
    1,
                    #   2
                      ]
shuffler.Delays    = [
    25.0001,
    # 15.0
    ]   # ns, adjust per source
shuffler.CollectionsSiPad = [
    "SiPadHits",
    # "SiPadHits"
]
shuffler.MaxEventsPerSource        = 50
shuffler.OutputFile                = "shuffled.edm4hep.root"
shuffler.OutputCollectionSiPad     = "SiPadHitsMerged"

ApplicationMgr(
    EvtSel  = "NONE",
    EvtMax  = 1,       # Only 1 Gaudi event: execute() is a no-op, all work is in finalize()
    TopAlg  = [shuffler],
    ExtSvc  = []
)

from k4FWCore import ApplicationMgr, IOSvc
from Configurables import EDM4HEP2RNTuple

SIPAD_BITFIELD = "system:8,layer:8,slice:5,x:9,y:9"

iosvc = IOSvc()
iosvc.Input = ["tracks.edm4hep.root"]

converter = EDM4HEP2RNTuple("EDM4HEP2RNTuple")
converter.InputFile      = "tracks.edm4hep.root"
converter.OutputFile     = "ShipHits.root"
converter.NTupleNames    = ["SiPad"]
converter.Collections    = ["SiPadHitsWindowed"]
converter.BitFields      = [SIPAD_BITFIELD]
converter.SourceIDParams = ["SiPadSourceIDs"]
converter.DetectorIDs    = [1]
converter.TrackFile           = "tracks.edm4hep.root"
converter.TrackCollectionName = "ACTSTracks"
converter.MeasCollections = ["SiPadMeasurements"]
converter.MeasNtupleNames = ["SiPadMeas"]
converter.MeasBitFields   = [SIPAD_BITFIELD]

ApplicationMgr(
    EvtSel  = "NONE",
    EvtMax  = -1,
    TopAlg  = [converter],
    ExtSvc  = [iosvc]
)

from k4FWCore import ApplicationMgr, IOSvc
from Configurables import ACTSGeoSvc, SiPadMeasConverter, ACTSProtoTracker
from Gaudi.Configuration import DEBUG, INFO

SIPAD_BITFIELD  = "system:8,layer:8,slice:5,x:9,y:9"
SIPAD_CELL_SIZE = 5.53  # mm, pixel pitch (Ecal_CellSizeX)

iosvc = IOSvc()
iosvc.Input          = "timewindows.edm4hep.root"
iosvc.Output         = "tracks.edm4hep.root"
iosvc.outputCommands = ["keep *"]

geosvc = ACTSGeoSvc("ACTSGeoSvc")
geosvc.CompactFile = "../../simulation/geometry/SND_compact.xml"
geosvc.OutputLevel = INFO

sipad_conv = SiPadMeasConverter("SiPadMeasConverter")
sipad_conv.GeoSvc           = "ACTSGeoSvc"
sipad_conv.InputCollection  = "SiPadHitsWindowed"
sipad_conv.OutputCollection = "SiPadMeasurements"
sipad_conv.BitField         = SIPAD_BITFIELD
sipad_conv.PixelSizeX       = SIPAD_CELL_SIZE
sipad_conv.PixelSizeY       = SIPAD_CELL_SIZE
sipad_conv.OutputLevel      = DEBUG

proto = ACTSProtoTracker("ACTSProtoTracker")
proto.GeoSvc           = "ACTSGeoSvc"
proto.InputSiPad       = "SiPadMeasurements"
proto.OutputCollection = "ACTSTracks"
proto.AutoSeed         = True
proto.MaxSeeds         = 3
proto.HoughBinSize     = 5.0
proto.HoughHalfSize    = 200.0
proto.HoughMinVotes    = 3
proto.SeedCompatRadius = 8.0
proto.SeedMomentum     = 3.0
proto.MaxChi2PerNdf    = 10.0
proto.HoughMaxMultiplicity = 10.0
proto.Chi2CutOff       = 70.0
proto.NumMeasCutOff    = 1
proto.IsolationWindow  = 5.0
proto.IsolationMaxNeighbors = 2
proto.OutputLevel      = DEBUG

ApplicationMgr(
    EvtSel  = "NONE",
    EvtMax  = -1,
    TopAlg  = [sipad_conv, proto],
    ExtSvc  = [iosvc, geosvc]
)

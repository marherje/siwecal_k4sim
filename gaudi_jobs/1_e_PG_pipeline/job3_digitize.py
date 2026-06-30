from k4FWCore import ApplicationMgr, IOSvc
from Configurables import GeV2MIPConversion, BasicDigitizer, DetectorFlipper, ChannelMapper
import os

REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))

infile = os.environ.get("INPUT_FILE", "timewindows.edm4hep.root")

iosvc = IOSvc()
iosvc.Input  = infile
iosvc.Output = "digitized.edm4hep.root"

mip = GeV2MIPConversion("GeV2MIP_SiPad")
mip.InputCollection  = "SiPadHitsWindowed"
mip.OutputCollection = "SiPadHitsMIP"
mip.MIPValue = 0.0002

dig = BasicDigitizer("BasicDigitizer_SiPad")
dig.InputCollection  = "SiPadHitsMIP"
dig.OutputCollection = "SiPadHitsDigi"
dig.Threshold = 0.0
dig.DebugFrequency = 500

flip = DetectorFlipper("DetectorFlipper_SiPad")
flip.InputCollection  = "SiPadHitsDigi"
flip.OutputCollection = "SiPadHitsFlipped"
flip.ZPositions = [
      0.0,   # slab  0  (2.8 mm W)
    -11.0,   # slab  1  (4.2 mm W)
    -22.0,   # slab  2
    -33.0,   # slab  3
    -44.0,   # slab  4
    -55.0,   # slab  5
    -66.0,   # slab  6
    -77.0,   # slab  7
    -88.0,   # slab  8
    -99.0,   # slab  9
   -110.0,   # slab 10
   -132.0,   # slab 11
   -143.0,   # slab 12
   -154.0,   # slab 13
   -165.0,   # slab 14
]
flip.BitField = "system:8,layer:8,slice:5,x:9,y:9"
flip.DebugFrequency = 500

cmap = ChannelMapper("ChannelMapper_SiPad")
cmap.InputCollection     = ["SiPadHitsFlipped"]
cmap.OutputCollection    = ["SiPadHitsMapped"]
cmap.OutputMaskedFlags   = ["SiPadHitsMasked"]
cmap.PadMapFile          = os.path.join(
    REPO_ROOT, "masking_info/geometry/fev10_rotate_chip_channel_x_y_mapping.txt")
cmap.PadMapFileSlab12    = os.path.join(
    REPO_ROOT, "masking_info/geometry/fev11_cob_good_rotate_chip_channel_x_y_mapping.txt")
cmap.MIPCalibFile        = os.path.join(
    REPO_ROOT, "masking_info/calibration/dummy_mip_map_15_highgain.txt")
cmap.MaxMIPValue         = 100.0
cmap.PositionTolerance   = 4.0
cmap.BitFieldIn          = "system:8,layer:8,slice:5,x:9,y:9"
cmap.BitFieldOut         = "system:8,slab:8,chip:16,channel:8,sca:8"
cmap.DebugFrequency      = 500

ApplicationMgr(
    EvtSel  = "NONE",
    EvtMax  = -1,
    TopAlg  = [mip, dig, flip, cmap],
    ExtSvc  = [iosvc]
)

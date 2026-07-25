from k4FWCore import ApplicationMgr, IOSvc
from Configurables import GeV2MIPConversion, BasicDigitizer, DetectorFlipper, ChannelMapper
import os
import yaml

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
with open(os.path.join(REPO_ROOT, "mappings", "slab_z_positions.yml")) as _zf:
    # Single source of truth for the per-slab z [mm], shared with the event
    # viewer and the compact geometry; never hardcode the array here.
    flip.ZPositions = [float(_z) for _z in yaml.safe_load(_zf)["slab_z_mm"]]
flip.BitField = "system:8,layer:8,slice:5,x:9,y:9"
flip.DebugFrequency = 500

cmap = ChannelMapper("ChannelMapper_SiPad")
cmap.InputCollection     = ["SiPadHitsFlipped"]
cmap.OutputCollection    = ["SiPadHitsMapped"]
cmap.OutputMaskedFlags   = ["SiPadHitsMasked"]
cmap.PadMapFile          = os.path.join(
    REPO_ROOT, "mappings/fev10_rotate_chip_channel_x_y_mapping.txt")
cmap.PadMapFileSlab12    = os.path.join(
    REPO_ROOT, "mappings/fev11_cob_good_rotate_chip_channel_x_y_mapping.txt")
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

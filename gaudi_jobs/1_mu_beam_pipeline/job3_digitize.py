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
mip.InputCollection  = "SiPadHits"
mip.OutputCollection = "SiPadHitsMIP"
# --- Single MIP value (scalar mode) ---
#mip.MIPValue = 0.0002
# --- Per-layer mode: uncomment and set MIPValues from mip_extraction_pipeline output ---
mip.MIPValues = [0.00020, 0.00015, 0.00015, 0.00015, 0.00015, 0.00015, 0.00015, 0.00015, 0.00015, 0.00015, 0.00015, 0.00015, 0.00015, 0.00015, 0.00015]

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

# ChannelMapper: translate simulation CellIDs to real TB channel IDs and apply
# MIP-calibration masking.
#
# - Reads SiPadHitsFlipped (sim Cartesian CellIDs)
# - Outputs SiPadHitsMapped (TB-format CellIDs: slab/chip/channel/sca)
# - Outputs SiPadHitsMasked (int32 UserDataCollection: 0=ok, 1=masked)
#
# MIP calibration: resolved from the MuonCalib_gaudi tree, see below.
cmap = ChannelMapper("ChannelMapper_SiPad")
cmap.InputCollection     = ["SiPadHitsFlipped"]
cmap.OutputCollection    = ["SiPadHitsMapped"]
cmap.OutputMaskedFlags   = ["SiPadHitsMasked"]
cmap.PadMapFile          = os.path.join(
    REPO_ROOT, "mappings/fev10_rotate_chip_channel_x_y_mapping.txt")
cmap.PadMapFileSlab12    = os.path.join(
    REPO_ROOT, "mappings/fev11_cob_good_rotate_chip_channel_x_y_mapping.txt")
# Masking comes from the muon calibration tree: mips/<threshold>/MIP_*_<gain>.txt.
# th230 masks ~3.5% of the channels (four dead chips in slab 0, one each in slabs
# 6 and 13, the rest scattered).  Switch CalibThreshold to th210/th220 for another
# trigger threshold, or set cmap.MIPCalibFile to bypass the tree with an explicit
# file (e.g. masking_info/calibration/dummy_mip_map_15_highgain.txt, which masks
# nothing).
cmap.CalibDir            = os.path.join(
    REPO_ROOT, "masking_info/calibration/MuonCalib_gaudi")
cmap.CalibThreshold      = "th230"
cmap.CalibGain           = "highgain"
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

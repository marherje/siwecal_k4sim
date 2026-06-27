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
# --- Single MIP value (scalar mode) ---
#mip.MIPValue = 0.0002
# --- Per-layer mode: uncomment and set MIPValues from mip_extraction_pipeline output ---
mip.MIPValues = [0.00020, 0.00015, 0.00015, 0.00015, 0.00015, 0.00015, 0.00015, 0.00015, 0.00015, 0.00015, 0.00015, 0.00015, 0.00015, 0.00015, 0.00015]

# mip.NLayers   = SIPAD_NLAYERS
# mip.BitField  = SIPAD_BITFIELD

dig = BasicDigitizer("BasicDigitizer_SiPad")
dig.InputCollection  = "SiPadHitsMIP"
dig.OutputCollection = "SiPadHitsDigi"
dig.Threshold = 0.5
dig.DebugFrequency = 500

# DetectorFlipper: rewrite the z coordinate of each hit to a canonical
# per-layer z table.  This is needed for real TB data where the detector is
# physically flipped with respect to the simulation convention.
#
# For SIMULATION the defaults reproduce the DD4hep z values, so the output
# is physically identical to SiPadHitsDigi (useful to test the algorithm and
# to produce a single canonical collection for downstream analysis).
#
# For real TB data, override ZPositions with the desired z convention, e.g.:
#   flip.ZPositions = list(reversed([...]))
#
# The collection SiPadHitsFlipped is what analysis/sim_to_ecal_tree.py reads.
flip = DetectorFlipper("DetectorFlipper_SiPad")
flip.InputCollection  = "SiPadHitsDigi"
flip.OutputCollection = "SiPadHitsFlipped"
# Default: simulation z positions (layer 0 at front = most negative z).
# Override here for real TB (flipped) data.
flip.ZPositions = [
    -116.35,  # layer  0  (2.8 mm W)
     -99.75,  # layer  1  (4.2 mm W)
     -83.15,  # layer  2
     -66.55,  # layer  3
     -49.95,  # layer  4
     -33.35,  # layer  5
     -16.75,  # layer  6
      -0.15,  # layer  7
      16.45,  # layer  8  (5.6 mm W)
      33.05,  # layer  9
      49.65,  # layer 10
      77.25,  # layer 11
      93.85,  # layer 12
     110.45,  # layer 13
     126.98,  # layer 14
]
flip.BitField = "system:8,layer:8,slice:5,x:9,y:9"
flip.DebugFrequency = 500

# ChannelMapper: translate simulation CellIDs to real TB channel IDs and apply
# MIP-calibration masking.
#
# - Reads SiPadHitsFlipped (sim Cartesian CellIDs)
# - Outputs SiPadHitsMapped (TB-format CellIDs: slab/chip/channel/sca)
# - Outputs SiPadHitsMasked (int32 UserDataCollection: 0=ok, 1=masked)
#
# Default MIP calibration: dummy file (all mpv=20 → nothing masked).
# Switch to real calibration file for data-driven masking:
#   cmap.MIPCalibFile = os.path.join(
#       REPO_ROOT, "masking_info/calibration/MuonCalib_it2_corrected/mips/th220",
#       "MIP_pedestalsubmode1_TB2026CERN_run_000142_highgain.txt")
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

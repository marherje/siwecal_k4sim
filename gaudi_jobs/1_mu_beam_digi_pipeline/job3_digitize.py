from k4FWCore import ApplicationMgr, IOSvc
from Configurables import GeV2MIPConversion, RealDigitizer
import os

SIPAD_BITFIELD = "system:8,layer:8,slice:5,x:9,y:9"
SIPAD_NLAYERS  = 15

infile = os.environ.get("INPUT_FILE", "timewindows.edm4hep.root")

iosvc = IOSvc()
iosvc.Input  = infile
iosvc.Output = "digitized.edm4hep.root"

mip = GeV2MIPConversion("GeV2MIP_SiPad")
mip.InputCollection  = "SiPadHitsWindowed"
mip.OutputCollection = "SiPadHitsMIP"
# --- Single MIP value (scalar mode) ---
mip.MIPValue = 0.0002
# --- Per-layer mode: uncomment and set MIPValues from mip_extraction_pipeline output ---
# mip.MIPValues = [0.0002, 0.0002, 0.0002, 0.0002, 0.0002,
#                  0.0002, 0.0002, 0.0002, 0.0002, 0.0002,
#                  0.0002, 0.0002, 0.0002, 0.0002, 0.0002]
# mip.NLayers   = SIPAD_NLAYERS
# mip.BitField  = SIPAD_BITFIELD

dig = RealDigitizer("RealDigitizer_SiPad")
dig.InputCollection    = "SiPadHitsMIP"
dig.OutputCollection   = "SiPadHitsDigi"
dig.Threshold = 0.0
dig.DigitizationMode   = "simple"   # "simple" = MIP cut (same as BasicDigitizer)
                                    # "real"   = realistic digitization (placeholder)
dig.DebugFrequency     = 500

ApplicationMgr(
    EvtSel  = "NONE",
    EvtMax  = -1,
    TopAlg  = [mip, dig],
    ExtSvc  = [iosvc]
)

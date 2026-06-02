from k4FWCore import ApplicationMgr, IOSvc
from Configurables import GeV2MIPConversion, BasicDigitizer
from Gaudi.Configuration import DEBUG
import os

infile = os.environ.get("INPUT_FILE", "")
if not infile:
    print("Use: INPUT_FILE=<input_file> k4run job3_digitize.py")
    import sys; sys.exit(1)

outfile = "digitized.edm4hep.root"

iosvc = IOSvc()
iosvc.Input = infile
iosvc.Output = outfile

mip = GeV2MIPConversion("GeV2MIP_SiPad")
mip.InputCollection  = "SiPadHitsWindowed"
mip.OutputCollection = "SiPadHitsMIP"
mip.MIPValue = 0.0002

dig = BasicDigitizer("BasicDigitizer_SiPad")
dig.InputCollection  = "SiPadHitsMIP"
dig.OutputCollection = "SiPadHitsDigi"
dig.Threshold = 0.5
dig.DebugFrequency = 500

ApplicationMgr(
    EvtSel  = "NONE",
    EvtMax  = -1,
    TopAlg  = [mip, dig],
    ExtSvc  = [iosvc]
)

"""
Step 0 of the digitisation: the per-layer MIP scale of the simulated stack.

MIPExtractor histograms the raw Geant4 energy per hit (GeV) layer by layer and
fits the MIP peak, so that GeV2MIPConversion / RealDigitizer can turn deposits
into MIP units with the value each layer really has -- a 650 um sensor collects
1.3x the charge of a 500 um one, and only the simulation itself knows what the
geometry file says.

Environment:
    INPUT_FILE        muon simulation to fit (default: timewindows.edm4hep.root,
                      the shuffled/split file of the old three-step pipeline; a
                      plain ddsim output works just as well and is what
                      mip_extraction_pipeline.sh feeds it)
    INPUT_COLLECTION  hit collection in that file (default: SiPadHitsWindowed for
                      the old pipeline, SiPadHits for a ddsim file)
    MIP_VALUES_OUT    text output (default: mip_values.txt); the shell script turns
                      it into mappings/mip_values_sim.yml, which job3 reads
    MIP_FIT_MODE      1 = Gaussian, 2 = Landau, 3 = Landau x Gaussian (default 3)
    MAX_EVENTS        events to read (default: all)
"""
import os

from k4FWCore import ApplicationMgr, IOSvc
from Configurables import MIPExtractor
from Gaudi.Configuration import INFO

SIPAD_BITFIELD = "system:8,layer:8,slice:5,x:9,y:9"
SIPAD_NLAYERS  = 15

infile = os.environ.get("INPUT_FILE", "timewindows.edm4hep.root")
collection = os.environ.get("INPUT_COLLECTION",
                            "SiPadHitsWindowed" if infile.endswith("timewindows.edm4hep.root")
                            else "SiPadHits")

iosvc = IOSvc()
iosvc.Input = infile

extractor = MIPExtractor("MIPExtractor")
extractor.InputCollection = collection            # raw hits in GeV (before GeV2MIP)
extractor.BitField        = SIPAD_BITFIELD
extractor.NLayers         = SIPAD_NLAYERS
extractor.FitMode         = int(os.environ.get("MIP_FIT_MODE", "3"))
extractor.MinEntries      = 100                    # minimum hits per layer to attempt fit
extractor.OutputRootFile  = os.environ.get("MIP_HISTOS_OUT", "mip_extraction.root")
extractor.OutputTextFile  = os.environ.get("MIP_VALUES_OUT", "mip_values.txt")
extractor.OutputLevel     = INFO

ApplicationMgr(
    EvtSel  = "NONE",
    EvtMax  = int(os.environ.get("MAX_EVENTS", "-1")),
    TopAlg  = [extractor],
    ExtSvc  = [iosvc]
)

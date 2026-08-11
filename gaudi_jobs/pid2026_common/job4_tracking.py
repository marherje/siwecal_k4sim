"""ACTS tracking over the SiW-ECAL silicon layers — shared by every pipeline.

Replaces the seven near-identical copies that used to live in the
gaudi_jobs/1_*_pipeline/ directories. Everything that varies between pipelines
comes from environment variables, everything that describes the detector comes
from the compact XML via parse_geometry:

    INPUT_FILE        input edm4hep file            (default timewindows.edm4hep.root)
    OUTPUT_FILE       output edm4hep file           (default tracks.edm4hep.root)
    INPUT_COLLECTION  hit collection to track on    (default SiPadHitsWindowed)
    SEED_MOMENTUM     beam momentum [GeV]           (default 3.0)
    TRACKING_LOGLEVEL DEBUG or INFO                 (default INFO)

Which INPUT_COLLECTION to use:
  * pipelines with time windows -> SiPadHitsWindowed on timewindows.edm4hep.root
  * per-chunk production        -> SiPadHitsDigi     on digitized.edm4hep.root

Use the PRE-flip collection (SiPadHitsDigi, not SiPadHitsFlipped/Mapped):
DetectorFlipper rewrites the hit z into the test-beam frame, which no longer
matches the ACTS surfaces — those live in the simulation frame that ACTSGeoSvc
reads from the same compact XML.
"""

import os
import sys
from pathlib import Path

from k4FWCore import ApplicationMgr, IOSvc
from Configurables import ACTSGeoSvc, SiPadMeasConverter, ACTSProtoTracker
from Gaudi.Configuration import DEBUG, INFO

REPO_ROOT = Path(__file__).resolve().parents[2]
GEOMETRY_DIR = REPO_ROOT / "simulation" / "geometry"
sys.path.insert(0, str(GEOMETRY_DIR))
from parse_geometry import SiWEcalGeometry  # noqa: E402

COMPACT_FILE = GEOMETRY_DIR / "SND_compact.xml"
geo = SiWEcalGeometry(COMPACT_FILE)

infile    = os.environ.get("INPUT_FILE", "timewindows.edm4hep.root")
outfile   = os.environ.get("OUTPUT_FILE", "tracks.edm4hep.root")
inputcoll = os.environ.get("INPUT_COLLECTION", "SiPadHitsWindowed")
seed_p    = float(os.environ.get("SEED_MOMENTUM", "3.0"))
loglevel  = DEBUG if os.environ.get("TRACKING_LOGLEVEL", "INFO") == "DEBUG" else INFO

iosvc = IOSvc()
iosvc.Input          = infile
iosvc.Output         = outfile
iosvc.outputCommands = ["keep *"]

# Builds the ACTS TrackingGeometry from the DD4hep DetElement tree: one plane
# surface per layer, at the sensitive slice, carrying that layer's full material
# stack. Absolute path so the job runs from any working directory.
geosvc = ACTSGeoSvc("ACTSGeoSvc")
geosvc.CompactFile = str(COMPACT_FILE)
geosvc.OutputLevel = loglevel

sipad_conv = SiPadMeasConverter("SiPadMeasConverter")
sipad_conv.GeoSvc           = "ACTSGeoSvc"
sipad_conv.InputCollection  = inputcoll
sipad_conv.OutputCollection = "SiPadMeasurements"
# Straight from the <readout> block, so it cannot drift from the segmentation.
sipad_conv.BitField         = geo.bitfields["SiPadHits"]
sipad_conv.PixelSizeX       = geo.ecal_cell_size_x
sipad_conv.PixelSizeY       = geo.ecal_cell_size_y
sipad_conv.OutputLevel      = loglevel

proto = ACTSProtoTracker("ACTSProtoTracker")
proto.GeoSvc           = "ACTSGeoSvc"
proto.InputSiPad       = "SiPadMeasurements"
proto.OutputCollection = "ACTSTracks"

# ---- Seeding (Hough transform over the 2D pad hits) ------------------------
proto.AutoSeed         = True
proto.MaxSeeds         = 3
proto.HoughBinSize     = 5.0
# Transverse envelope from the geometry (was hardcoded to 200 mm, inherited
# from the SND setup; this detector is 180 mm across).
proto.HoughHalfSize    = geo.hough_half_size
proto.HoughMinVotes    = 3
proto.SeedCompatRadius = 8.0
proto.SeedMomentum     = seed_p
# Refinement/regression bin for the seed position; one pad pitch.
proto.SeedStripPitch   = geo.ecal_cell_size_x

# ---- Track finding --------------------------------------------------------
proto.Chi2CutOff       = 70.0
proto.NumMeasCutOff    = 1
proto.MaxChi2PerNdf    = 10.0

# ---- Shower rejection ----------------------------------------------------
# Both neighbour-counting windows MUST be expressed in pad pitches. Any window
# <= one pitch counts zero neighbours by construction (the nearest other pad on
# a plane is exactly one pitch away), which silently disabled both filters: the
# 5.0 mm inherited here came from the SND SiTarget, where the strips are 75 um
# and 5 mm spans dozens of channels. 1.5 pitches reaches the 4 orthogonal
# neighbours (1 pitch) and the 4 diagonal ones (1.41 pitches) — the 8 pads
# surrounding a hit — and nothing beyond.
NEIGHBOUR_WINDOW = 1.5 * geo.ecal_cell_size_x

# Hough-level: a peak fed by many crossings per layer is a shower, not a track.
proto.HoughMaxMultiplicity  = 10.0
proto.IsolationWindow       = NEIGHBOUR_WINDOW
proto.IsolationMaxNeighbors = 2
# Pool-level: drop a hit when more than half of the 8 pads around it are also
# lit, i.e. it sits in a shower core rather than on a MIP trail. This is what
# keeps an electron shower from dragging the incoming track off.
proto.HitPurgeWindow        = NEIGHBOUR_WINDOW
proto.HitPurgeMaxNeighbors  = 4

# ---- Track quality -------------------------------------------------------
proto.SeedCleaning            = True
proto.FinalRefit              = True
proto.DuplicateOverlapFraction = 0.7

proto.OutputLevel      = loglevel

ApplicationMgr(
    EvtSel  = "NONE",
    EvtMax  = -1,
    TopAlg  = [sipad_conv, proto],
    ExtSvc  = [iosvc, geosvc]
)

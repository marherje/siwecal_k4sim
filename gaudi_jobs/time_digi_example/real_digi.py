from k4FWCore import ApplicationMgr, IOSvc
from Configurables import RealDigitizer
import os
from pathlib import Path


def env_float(name, default):
    return float(os.environ.get(name, default))


def env_int(name, default):
    return int(os.environ.get(name, default))


default_input = Path(
    "/home/llr/ilc/shi/data/siwecal_k4sim/output/"
    "output_PG_gamma_10GeV_100evt.edm4hep.root"
)
infile = Path(os.environ.get("INPUT_FILE", default_input))
default_output = infile.with_name(
    infile.name.replace(".edm4hep.root", "_real_digitized.edm4hep.root")
)
outfile = Path(os.environ.get("OUTPUT_FILE", default_output))

iosvc = IOSvc()
iosvc.Input = str(infile)
iosvc.Output = str(outfile)

dig = RealDigitizer("RealDigitizer_SiPad")
dig.InputCollection = os.environ.get("INPUT_COLLECTION", "SiPadHits")
dig.OutputCollection = os.environ.get("OUTPUT_COLLECTION", "SiPadHitsDigi")
dig.DigitizedEnergyCollection = os.environ.get(
    "DIGITIZED_ENERGY_COLLECTION", "SiPadHitsDigiDigitizedEnergy"
)
dig.DigitizedTimeCollection = os.environ.get(
    "DIGITIZED_TIME_COLLECTION", "SiPadHitsDigiDigitizedTime"
)
dig.InputEnergyUnit = "GeV"
dig.MIPValue = env_float("MIP_VALUE_GEV", 0.0002)
dig.Threshold = env_float("MIP_THRESHOLD", 0.5)
dig.DigitizationMode = "real"
dig.DebugFrequency = env_int("DEBUG_FREQUENCY", 1)

dig.DelayNs = env_float("DELAY_NS", 160.0)
dig.TauFastNs = env_float("TAU_FAST_NS", 30.0)
dig.TauSlowNs = env_float("TAU_SLOW_NS", 180.0)
dig.OrderFast = env_int("ORDER_FAST", 2)
dig.OrderSlow = env_int("ORDER_SLOW", 2)
dig.FastWindowNs = env_float("FAST_WINDOW_NS", 200.0)
dig.SlowWindowNs = env_float("SLOW_WINDOW_NS", 500.0)
dig.FastNoiseMIP = env_float("FAST_NOISE_MIP", 1.0 / 30.0)
dig.SlowNoiseMIP = env_float("SLOW_NOISE_MIP", 1.0 / 12.0)
dig.PeakSearchBins = env_int("PEAK_SEARCH_BINS", 64)
dig.RefineIterations = env_int("REFINE_ITERATIONS", 48)
dig.TriggerSearchBins = env_int("TRIGGER_SEARCH_BINS", 64)
dig.RandomSeed = env_int("RANDOM_SEED", 5489)

ApplicationMgr(
    EvtSel="NONE",
    EvtMax=-1,
    TopAlg=[dig],
    ExtSvc=[iosvc],
)

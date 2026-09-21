"""Job 3 — digitisation, z flip and channel mapping.

Shared by the condor chunk jobs (simulation/run_script/generic_condor_beam_chunk.sh,
reprocess_chunk.sh) and by the local 1_* pipelines.

CALIB_THRESHOLD (th210 / th220 / th230, default th230) selects the threshold set:
it picks the masking table AND the trigger model of the real chain, both measured
from the test beam.  DIGI_MODE selects which digitisation chain(s) run:

  simple (default)  SiPadHits -> GeV2MIPConversion -> BasicDigitizer
                              -> DetectorFlipper   -> ChannelMapper
                              => SiPadHitsMapped + SiPadHitsMasked
                    This is what production has always run; leaving DIGI_MODE
                    unset changes nothing.

  real              SiPadHits -> RealDigitizer(mode='real')
                              -> DetectorFlipper -> ChannelMapper
                              -> AdcDigitizer (unless ADC_MODEL=0)
                              => SiPadHitsRealAdc + SiPadHitsRealAdcHigh/Low
                                 (+ SiPadHitsRealMapped/Masked, ...Time, ...Fast)

  both              Both chains on the same input, into the same output file, so
                    the two can be compared hit by hit without re-simulating.

The real chain reads SiPadHits in GeV *before* GeV2MIPConversion on purpose:
GeV2MIPConversion (like BasicDigitizer, DetectorFlipper and ChannelMapper)
creates fresh hits without copying CaloHitContributions, and the cell shaping
needs them.  RealDigitizer normalises to MIP itself via MIPValues, so
GeV2MIPConversion would be redundant on that branch anyway.
"""

from k4FWCore import ApplicationMgr, IOSvc
from Configurables import (GeV2MIPConversion, BasicDigitizer, DetectorFlipper,
                           ChannelMapper, RealDigitizer, AdcDigitizer)
import os
import sys
import yaml

REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))

SIPAD_BITFIELD_SIM = "system:8,layer:8,slice:5,x:9,y:9"
SIPAD_BITFIELD_TB  = "system:8,slab:8,chip:16,channel:8,sca:8"
SIPAD_NLAYERS      = 15

# Per-layer MIP calibration [GeV/MIP]: the peak of the raw energy per hit of a
# simulated muon sample, layer by layer, as gaudi_jobs/mip_extraction_pipeline/
# writes it to mappings/mip_values_sim.yml.  Step 0 of the digitisation, never a
# list typed here: the value depends on the geometry (a 650 um sensor is 1.34x a
# 500 um one, and the old hand-pasted list had that factor on the wrong layer).
# Defined once, because GeV2MIPConversion (simple chain) and RealDigitizer (real
# chain) must agree or the two chains are not comparable.
MIP_VALUES_FILE = os.environ.get("MIP_VALUES_FILE",
                                 os.path.join(REPO_ROOT, "mappings", "mip_values_sim.yml"))


def _load_mip_values(path):
    if not os.path.isfile(path):
        raise SystemExit(
            f"MIP values file not found: {path}\n"
            "Run gaudi_jobs/mip_extraction_pipeline/mip_extraction_pipeline.sh <muon ddsim> "
            "first (step 0), or point MIP_VALUES_FILE at its output.")
    with open(path) as fh:
        doc = yaml.safe_load(fh)
    values = [float(v) for v in doc.get("mip_gev", [])]
    if len(values) != SIPAD_NLAYERS or any(v <= 0 for v in values):
        raise SystemExit(f"{path}: expected {SIPAD_NLAYERS} positive per-layer MIP values, got {values}")
    return values


MIP_VALUES = _load_mip_values(MIP_VALUES_FILE)

# Two MIP units, on purpose.  The SIMPLE chain wants "one crossing of THIS layer
# = 1 MIP", i.e. the per-layer values above, because its hit_energy is the
# analysis quantity and the data's reconstruction normalises every layer by its
# own table.  The REAL chain wants a CHARGE unit: the ADC is proportional to the
# collected charge, so a 650 um crossing must give 1.34x the ADC of a 500 um one,
# and the reconstruction table (per layer) is what turns that back into the
# data's MIP.  Its unit is therefore the crossing of a reference (500 um) layer,
# the same for every layer, and the trigger threshold (threshold_adc /
# adc_per_mip) lives in that unit too.  MIP_VALUES_CHARGE is set below once the
# slab description says which layers are the reference ones.

# --------------------------------------------------------------------------- #
# Threshold set.  ONE input drives everything that depends on it: the dead-channel
# masking (ChannelMapper reads mips/<th>/) and the trigger model of the real
# chain (threshold, spread and plateau efficiency, measured per th from the test
# beam's own hitbit_high -- see analysis/trigger_turnon.py).
#
# Before this, RealDigitizer applied a hardcoded 0.5 MIP step with 100%
# efficiency whatever the sample was meant to represent, while the measured
# discriminator sits at 16.3 / 19.7 / 25.4 ADC for th210 / th220 / th230 and
# plateaus at 85-96%.
# --------------------------------------------------------------------------- #
CALIB_THRESHOLD = os.environ.get("CALIB_THRESHOLD", "th230")

with open(os.path.join(REPO_ROOT, "mappings", "digi_calibration.yml")) as _cf:
    _CALIB = yaml.safe_load(_cf)["thresholds"]
if CALIB_THRESHOLD not in _CALIB:
    raise SystemExit(
        f"CALIB_THRESHOLD='{CALIB_THRESHOLD}' has no entry in "
        f"mappings/digi_calibration.yml (have: {', '.join(sorted(_CALIB))}). "
        f"Measure its trigger with analysis/trigger_turnon.py and its ADC scale "
        f"with analysis/fit_adc_scale.py.")
TH = _CALIB[CALIB_THRESHOLD]

# The ADC model closes the chain: MIP -> ADC with THIS threshold set's own fitted
# scale, then back to MIP through the SAME calibration table the data run was
# reconstructed with.  Without it the simulation divides by its own clean MIP
# scale while the data divides by a measured (and, at th230, biased) one, so the
# two energies are not the same quantity.  ADC_MODEL=0 goes back to the raw
# digitised MIP.
ADC_MODEL = os.environ.get("ADC_MODEL", "1") not in ("0", "no", "false")
# Only the channel-to-channel VARIATION of the gain comes from a table; the
# absolute scale is TH["adc_per_mip"], fitted against this threshold's own data.
GAIN_SHAPE_THRESHOLD = os.environ.get("GAIN_SHAPE_THRESHOLD", "th210")

HIT_SELECTION = os.environ.get("HIT_SELECTION", "cell")
if HIT_SELECTION not in ("cell", "chip"):
    raise SystemExit(f"HIT_SELECTION='{HIT_SELECTION}' is not one of 'cell', 'chip'.")
ADC_HIT_THRESHOLD = float(os.environ.get("ADC_HIT_THRESHOLD", "30"))
DIGI_MODE = os.environ.get("DIGI_MODE", "simple")
if DIGI_MODE not in ("simple", "real", "both"):
    raise SystemExit(
        f"DIGI_MODE='{DIGI_MODE}' is not one of 'simple', 'real', 'both'.")

infile = os.environ.get("INPUT_FILE", "timewindows.edm4hep.root")

iosvc = IOSvc()
iosvc.Input  = infile
iosvc.Output = "digitized.edm4hep.root"
# Explicit so the parallel UserDataCollections written by RealDigitizer
# (digitised energy, trigger time) do not depend on an IOSvc default.
iosvc.outputCommands = ["keep *"]

# Single source of truth for the per-slab z [mm], shared with the event viewer
# and the compact geometry; never hardcode the array here.
with open(os.path.join(REPO_ROOT, "mappings", "slab_z_positions.yml")) as _zf:
    SLAB_Z = [float(_z) for _z in yaml.safe_load(_zf)["slab_z_mm"]]

COMPACT_FILE = os.path.join(REPO_ROOT, "simulation", "geometry", "SND_compact.xml")
CALIB_DIR    = os.path.join(REPO_ROOT, "masking_info/calibration/MuonCalib_gaudi")

# The per-slab hardware description (technology, sensor thickness, pad map,
# threshold DAC) -- the same YAML the test-beam repo carries.  Which slab is the
# chip-on-board and needs its own pad map is read from here, not typed.
sys.path.insert(0, REPO_ROOT)
from analysis.slab_description import load_slab_description  # noqa: E402
SLABS = load_slab_description(os.path.join(REPO_ROOT, "mappings", "slab_z_positions.yml"))
PAD_MAP = SLABS.default_pad_map()
PAD_MAP_OVERRIDES = [f"{slab}:{path}" for slab, path in sorted(SLABS.pad_map_overrides().items())]
_ref = [v for v, t in zip(MIP_VALUES, SLABS.sensor_thickness_um) if t == 500]
MIP_CHARGE_GEV = sum(_ref) / len(_ref) if _ref else sum(MIP_VALUES) / len(MIP_VALUES)
MIP_VALUES_CHARGE = [MIP_CHARGE_GEV] * SIPAD_NLAYERS


def make_flipper(name, in_coll, out_coll):
    """DetectorFlipper: hit z -> test-beam frame."""
    flip = DetectorFlipper(name)
    flip.InputCollection  = in_coll
    flip.OutputCollection = out_coll
    flip.ZPositions = SLAB_Z
    # Lets DetectorFlipper cross-check the YAML's layer spacing against the geometry
    # and warn if the two have drifted apart -- the check that would have caught the
    # table being left at the old 11 mm pitch.
    flip.CompactFile = COMPACT_FILE
    flip.BitField = SIPAD_BITFIELD_SIM
    flip.DebugFrequency = 500
    return flip


def make_mapper(name, in_coll, out_coll, masked_coll):
    """ChannelMapper: sim CellIDs -> TB channel IDs + dead-channel masking.

    Masking comes from the muon calibration tree: mips/<threshold>/MIP_*_<gain>.txt.
    th230 masks ~3.5% of the channels (four dead chips in slab 0, one each in slabs
    6 and 13, the rest scattered).  Switch CalibThreshold to th210/th220 for another
    trigger threshold, or set MIPCalibFile to bypass the tree with an explicit file
    (e.g. masking_info/calibration/dummy_mip_map_15_highgain.txt, which masks
    nothing).
    """
    cmap = ChannelMapper(name)
    cmap.InputCollection   = [in_coll]
    cmap.OutputCollection  = [out_coll]
    cmap.OutputMaskedFlags = [masked_coll]
    cmap.PadMapFile        = PAD_MAP
    cmap.PadMapSlabOverrides = PAD_MAP_OVERRIDES
    cmap.CalibDir          = CALIB_DIR
    cmap.CalibThreshold    = CALIB_THRESHOLD
    cmap.CalibGain         = "highgain"
    cmap.MaxMIPValue       = 100.0
    cmap.PositionTolerance = 4.0
    cmap.BitFieldIn        = SIPAD_BITFIELD_SIM
    cmap.BitFieldOut       = SIPAD_BITFIELD_TB
    cmap.DebugFrequency    = 500
    return cmap


top_alg = []

# --------------------------------------------------------------------------- #
# Simple chain — the production default
# --------------------------------------------------------------------------- #
if DIGI_MODE in ("simple", "both"):
    mip = GeV2MIPConversion("GeV2MIP_SiPad")
    mip.InputCollection  = "SiPadHits"
    mip.OutputCollection = "SiPadHitsMIP"
    mip.MIPValues = MIP_VALUES
    mip.NLayers   = SIPAD_NLAYERS
    mip.BitField  = SIPAD_BITFIELD_SIM

    dig = BasicDigitizer("BasicDigitizer_SiPad")
    dig.InputCollection  = "SiPadHitsMIP"
    dig.OutputCollection = "SiPadHitsDigi"      # <-- tracking input (pre-flip)
    dig.Threshold = 0.0
    dig.DebugFrequency = 500

    top_alg += [mip, dig,
                make_flipper("DetectorFlipper_SiPad",
                             "SiPadHitsDigi", "SiPadHitsFlipped"),
                make_mapper("ChannelMapper_SiPad",
                            "SiPadHitsFlipped", "SiPadHitsMapped",
                            "SiPadHitsMasked")]

# --------------------------------------------------------------------------- #
# Real chain — contribution-level CR-RC cell shaping
# --------------------------------------------------------------------------- #
if DIGI_MODE in ("real", "both"):
    real = RealDigitizer("RealDigitizer_SiPad")
    real.InputCollection  = "SiPadHits"          # GeV, with contributions
    real.OutputCollection = "SiPadHitsRealDigi"  # <-- tracking input (pre-flip)
    real.DigitizedEnergyCollection = "SiPadHitsRealDigitizedEnergy"
    real.DigitizedTimeCollection   = "SiPadHitsRealDigitizedTime"
    real.DigitizedFastCollection   = "SiPadHitsRealDigitizedFast"
    real.DigitizedTriggerCollection = "SiPadHitsRealDigitizedTrigger"
    # HIT_SELECTION=chip mirrors EcalEventBuilder's HitSelection='adc': every
    # sampled cell is written with its trigger flag and AdcDigitizer keeps, on
    # chips where a cell fired, the cells that fired or exceed ADC_HIT_THRESHOLD.
    # 'cell' (default) is the historical behaviour, one hit per fired cell.
    real.HitSelection = HIT_SELECTION
    real.DigitizationMode = "real"
    real.InputEnergyUnit  = "GeV"
    real.HitEnergyContent = "digitized"          # shaped slow sample [MIP] -> hit energy
    real.MIPValues = MIP_VALUES_CHARGE   # charge unit: see MIP_VALUES above
    real.NLayers   = SIPAD_NLAYERS
    real.BitField  = SIPAD_BITFIELD_SIM
    # The trigger model of this threshold set.  Threshold is the discriminator
    # level, FastNoiseMIP the spread that gives the turn-on its width (noise plus
    # channel-to-channel dispersion, indistinguishable here), and
    # TriggerEfficiency the plateau above it.  REAL_THRESHOLD_MIP still overrides
    # the level by hand for a scan.
    # The discriminator level is measured in ADC; what it is worth in MIP follows
    # from this set's own fitted scale.  Both numbers come from the same YAML
    # entry, so the threshold moves with the scale whenever the fit is redone.
    _adc_per_mip = float(os.environ.get("ADC_PER_MIP_OVERRIDE",
                                        TH["adc_per_mip"]))
    _threshold_mip = TH["threshold_adc"] / _adc_per_mip
    real.Threshold = float(os.environ.get("REAL_THRESHOLD_MIP", _threshold_mip))
    real.FastNoiseMIP = float(os.environ.get("REAL_FAST_NOISE_MIP",
                                             TH["sigma_mip"]))
    real.TriggerEfficiency = float(os.environ.get("REAL_TRIGGER_EFF",
                                                  TH["efficiency"]))
    # Slabs that ran at a discriminator of their own (slab 12, the chip-on-board,
    # at DAC 243 against the set's 215-230): digi_calibration.yml carries their
    # measured turn-on under `slab_overrides`, in ADC, and it is worth
    # threshold_adc / adc_per_mip of the SET -- their gain is normal, only the
    # discriminator differs.  The env overrides above apply to every layer alike.
    _overrides = TH.get("slab_overrides") or {}
    if _overrides and "REAL_THRESHOLD_MIP" not in os.environ:
        _thr = [real.Threshold] * SIPAD_NLAYERS
        _sig = [real.FastNoiseMIP] * SIPAD_NLAYERS
        _eff = [real.TriggerEfficiency] * SIPAD_NLAYERS
        for _slab, _o in _overrides.items():
            _slab = int(_slab)
            if "threshold_adc" in _o:
                _thr[_slab] = float(_o["threshold_adc"]) / _adc_per_mip
            if "sigma_adc" in _o:
                _sig[_slab] = float(_o["sigma_adc"]) / _adc_per_mip
            if "efficiency" in _o:
                _eff[_slab] = float(_o["efficiency"])
        real.ThresholdPerLayer = _thr
        real.FastNoiseMIPPerLayer = _sig
        real.TriggerEfficiencyPerLayer = _eff
    real.RandomSeed = int(os.environ.get("REAL_RANDOM_SEED", 5489))
    real.DebugFrequency = 500

    top_alg += [real,
                make_flipper("DetectorFlipper_SiPadReal",
                             "SiPadHitsRealDigi", "SiPadHitsRealFlipped"),
                make_mapper("ChannelMapper_SiPadReal",
                            "SiPadHitsRealFlipped", "SiPadHitsRealMapped",
                            "SiPadHitsRealMasked")]

    if ADC_MODEL:
        # After ChannelMapper, which is what knows the test-beam channel and so
        # which pedestal and MPV apply.  Strictly 1:1, so the parallel masking /
        # time / fast-peak collections stay index-aligned with the hits.
        adc = AdcDigitizer("AdcDigitizer_SiPad")
        adc.InputCollection   = "SiPadHitsRealMapped"
        adc.OutputCollection  = "SiPadHitsRealAdc"
        adc.AdcHighCollection = "SiPadHitsRealAdcHigh"
        adc.AdcLowCollection  = "SiPadHitsRealAdcLow"
        adc.CalibDir        = CALIB_DIR
        adc.CalibThreshold      = CALIB_THRESHOLD        # the table that RECONSTRUCTS
        # ADC_PER_MIP_OVERRIDE is how analysis/fit_adc_scale.py turns the scale
        # while it iterates, without rewriting the YAML on every pass.
        adc.AdcPerMip           = float(os.environ.get("ADC_PER_MIP_OVERRIDE",
                                                       TH["adc_per_mip"]))
        adc.GainShapeThreshold  = GAIN_SHAPE_THRESHOLD   # shape only, never the scale
        # The dynamic range of THIS set, fitted from its electron data
        # (analysis/compare_gains.py --write): the roll-over's ceiling and
        # sharpness, and the LG<->HG anchor line the reconstruction inverts.
        adc.AdcHighMax          = float(TH["adc_high_max"])
        adc.SaturationOrder     = float(TH.get("saturation_order", 0.0))
        adc.SaturationAdc       = float(TH["saturation_adc"])
        adc.GainRatio           = float(TH.get("gain_ratio", 0.0961))
        adc.GainIntercept       = float(TH.get("gain_intercept", 1.60))
        # The low gain's own scatter about that line (compare_gains.py --write):
        # pedestal width plus the gain-ratio dispersion; widens the energy of
        # the hits reconstructed from the low gain, moves nothing on average.
        adc.LowGainNoiseAdc     = float(TH.get("lg_noise_adc", 0.0))
        adc.LowGainSpread       = float(TH.get("lg_gain_spread", 0.0))
        adc.CalibGain           = "highgain"
        # per-slab: a slab's table offset (slab 12's truncated MIP fit, layer 14's
        # thicker sensor that the geometry already carries) is not gain.
        adc.GainShapeNormalisation = os.environ.get("GAIN_SHAPE_NORM", "per-slab")
        if HIT_SELECTION == "chip":
            adc.TriggerCollection = "SiPadHitsRealDigitizedTrigger"
            adc.AdcHitThreshold   = ADC_HIT_THRESHOLD
        adc.KeptCollection      = "SiPadHitsRealAdcKept"
        adc.BitField        = SIPAD_BITFIELD_TB
        adc.DebugFrequency  = 500
        top_alg += [adc]

ApplicationMgr(
    EvtSel  = "NONE",
    EvtMax  = -1,
    TopAlg  = top_alg,
    ExtSvc  = [iosvc]
)

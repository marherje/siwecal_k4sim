Digitiser table tree for the final_v5 simulation (DIGI_CALIB_DIR of gaudi_jobs/pid2026_common/job3_digitize.py).

Every threshold slot holds the th210 MIP table re-derived on the fixed-SCA-pairing chunks
(siwecal-tb2026 calibration/MuonCalib_gaudi_fixed/mips/th210, eudaq 152-155, 158-160), because the data campaign
it is compared with (Reconstructed_final) calibrates AND masks every run with th210. So in the simulation:
  - ChannelMapper masking = th210 table for every set (same dead/uncalibrated channels as the data);
  - AdcDigitizer gain shape and reconstruction table = th210 (RECO_TABLE_THRESHOLD=th210);
  - the trigger model (threshold, spread, plateau) still comes from mappings/digi_calibration.yml for the set
    (CALIB_THRESHOLD), which is what differs between th210/th220/th230 runs.
The th230/th220 tables of the fixed calibration are not used here (th230 fits only 3452 channels: as a mask it
would switch off half the detector).

pedestals/th{210,220,230}: the fixed th210 pedestal table (same reason: the data use it for every run). With
PED_NOISE=1 the digitiser adds each channel's electronic noise from its widths (mean over the channel's fitted
SCAs, median 1.37 ADC high gain) and RealDigitizer's scalar SlowNoiseMIP is off.

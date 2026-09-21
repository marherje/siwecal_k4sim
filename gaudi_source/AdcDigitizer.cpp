#include "Gaudi/Algorithm.h"
#include "GaudiKernel/MsgStream.h"
#include "k4FWCore/DataHandle.h"
#include "CalibTables.hh"
#include "edm4hep/SimCalorimeterHitCollection.h"
#include "podio/UserDataCollection.h"
#include "DDSegmentation/BitFieldCoder.h"

#include <atomic>
#include <algorithm>
#include <cmath>
#include <iomanip>
#include <vector>
#include <cstdint>
#include <map>
#include <utility>
#include <filesystem>
#include <memory>
#include <random>
#include <string>

// AdcDigitizer — turns the digitised MIP amplitude into an ADC reading, and then
// reconstructs it back to MIP the way the test beam's reconstruction does.
//
// Why this exists
// ---------------
// Up to here the simulated chain works in MIP from end to end, so it never has
// the two things the real readout has: a finite ADC with a pedestal and a
// saturating high gain, and -- more important -- a CALIBRATION standing between
// the charge and the energy the analysis sees.
//
// The test beam's energy is `(adc_high - pedestal) / MPV`, where MPV comes from a
// muon run.  That MPV is not perfect: at the highest threshold the trigger eats
// into the MIP peak and the fitted MPV comes out high (measured: th230 sits 1.47x
// above th220 channel by channel, while the pedestals and the LG->HG anchor say
// the ADC-per-charge scale did NOT change).  A simulation that divides by its own
// clean MIP scale is therefore NOT comparable with data that divides by a biased
// one, however good its physics is.
//
// So this algorithm does what the detector plus the reconstruction do, in order:
//
//   1. MIP -> ADC with `AdcPerMip`, the ADC a channel gives per MIP.  This number
//      is MEASURED PER THRESHOLD SET, by anchoring the simulation to test-beam
//      data of that same set at a reference energy (analysis/fit_adc_scale.py),
//      and it is the only absolute scale in the chain.  It used to be lifted
//      wholesale from one threshold set's MIP table on the argument that that
//      table was the least biased -- an inference about which table is wrong,
//      doing real work in the answer.  Now nothing is inherited between sets.
//
//      Only the channel-to-channel VARIATION comes from a table
//      (`GainShapeThreshold`), normalised to its own median, because that is the
//      part the three sets agree on to ~11%; the scale is data.
//   2. + pedestal, round to integer ADC, clamp the high gain at `AdcHighMax`
//      (it saturates), and carry the low gain through the LG<->HG anchor line,
//      which stays linear over the whole range.
//   3. ADC -> MIP with the RECONSTRUCTION table (`CalibThreshold`), including the
//      low-gain recovery above `SaturationAdc` -- a port of
//      EventBuilder::buildHit, so the simulation inherits exactly the same
//      calibration, and the same bias, as the data it is compared with.
//
// A channel missing from either table gets energy 0, which is what
// EventBuilder::buildHit does with a masked channel.  The algorithm is strictly
// 1:1: every input hit produces one output hit, because the parallel
// UserDataCollections written upstream (masking flags, trigger time, fast peak)
// are tied to the hits by index alone.
class AdcDigitizer : public Gaudi::Algorithm {
public:
  AdcDigitizer(const std::string& name, ISvcLocator* svcLoc)
      : Gaudi::Algorithm(name, svcLoc) {}

  StatusCode initialize() override {
    try {
      m_inputHandle = std::make_unique<k4FWCore::DataHandle<edm4hep::SimCalorimeterHitCollection>>(
          m_inputName.value(), Gaudi::DataHandle::Reader, this);
      m_outputHandle = std::make_unique<k4FWCore::DataHandle<edm4hep::SimCalorimeterHitCollection>>(
          m_outputName.value(), Gaudi::DataHandle::Writer, this);
      m_adcHighHandle = std::make_unique<k4FWCore::DataHandle<podio::UserDataCollection<float>>>(
          m_adcHighName.value(), Gaudi::DataHandle::Writer, this);
      m_adcLowHandle = std::make_unique<k4FWCore::DataHandle<podio::UserDataCollection<float>>>(
          m_adcLowName.value(), Gaudi::DataHandle::Writer, this);
      m_keptHandle = std::make_unique<k4FWCore::DataHandle<podio::UserDataCollection<int>>>(
          m_keptName.value(), Gaudi::DataHandle::Writer, this);
      if (!m_triggerName.value().empty()) {
        m_triggerHandle = std::make_unique<k4FWCore::DataHandle<podio::UserDataCollection<int>>>(
            m_triggerName.value(), Gaudi::DataHandle::Reader, this);
        info() << "[AdcDigitizer] chip-level hit selection: a chip counts when one of its "
                  "cells fired (" << m_triggerName.value() << "); on such a chip a cell is a "
                  "hit if it fired or its high gain exceeds " << m_adcHitThreshold.value()
               << " ADC over the pedestal -- EcalEventBuilder HitSelection='adc'" << endmsg;
      }

      if (!loadTables()) return StatusCode::FAILURE;

      m_decoder = std::make_unique<dd4hep::DDSegmentation::BitFieldCoder>(m_bitField.value());

      info() << "[AdcDigitizer] AdcPerMip=" << m_adcPerMip.value()
             << " ADC/MIP (shape from " << m_gainShapeThreshold.value() << ", "
             << m_gainTable.nCalibrated << " channels), reconstruction table "
             << m_calibThreshold.value() << " (" << m_recoTable.nCalibrated
             << " channels), pedestals " << m_pedestals.nCalibrated << " entries"
             << endmsg;
      info() << "[AdcDigitizer] AdcHighMax=" << m_adcHighMax.value()
             << "  SaturationOrder=" << m_saturationOrder.value()
             << "  SaturationAdc=" << m_saturationAdc.value()
             << "  GainRatio=" << m_gainRatio.value()
             << "  GainIntercept=" << m_gainIntercept.value()
             << "  LowGainNoiseAdc=" << m_lowGainNoiseAdc.value()
             << "  LowGainSpread=" << m_lowGainSpread.value()
             << "  Quantise=" << (m_quantise.value() ? "yes" : "no") << endmsg;
      return Gaudi::Algorithm::initialize();
    } catch (const std::exception& e) {
      error() << "[AdcDigitizer] Exception in initialize(): " << e.what() << endmsg;
      return StatusCode::FAILURE;
    }
  }

  StatusCode execute(const EventContext&) const override {
    try {
      const auto* input = m_inputHandle->get();
      auto* output = m_outputHandle->createAndPut();
      auto* adcHighOut = m_adcHighHandle->createAndPut();
      auto* adcLowOut = m_adcLowHandle->createAndPut();
      auto* keptOut = m_keptHandle->createAndPut();

      // The per-cell trigger flags, if the chip-level selection is on.  With
      // them, which chips fired is known before the loop: a chip fired when at
      // least one of its cells did.  Index alignment is the only tie to the hits,
      // as for every other parallel collection of the chain.
      const podio::UserDataCollection<int>* trigger = nullptr;
      std::map<std::pair<int, int>, bool> chipFired;
      if (m_triggerHandle) {
        trigger = m_triggerHandle->get();
        if (trigger->size() != input->size()) {
          error() << "[AdcDigitizer] trigger collection has " << trigger->size()
                  << " entries for " << input->size() << " hits" << endmsg;
          return StatusCode::FAILURE;
        }
        std::size_t i = 0;
        for (const auto& hit : *input) {
          const std::uint64_t cellID = hit.getCellID();
          const auto key = std::make_pair(static_cast<int>(m_decoder->get(cellID, "slab")),
                                          static_cast<int>(m_decoder->get(cellID, "chip")));
          if ((*trigger)[i] > 0) chipFired[key] = true;
          ++i;
        }
      }

      const long long evtNum = m_eventCount.fetch_add(1);
      const bool doPrint = (evtNum % m_debugFreq == 0);
      int nUncalibrated = 0, nSaturated = 0, nDropped = 0, nRecovered = 0;

      std::size_t hitIndex = 0;
      for (const auto& hit : *input) {
        const std::uint64_t cellID = hit.getCellID();
        const int slab = static_cast<int>(m_decoder->get(cellID, "slab"));
        const int chip = static_cast<int>(m_decoder->get(cellID, "chip"));
        const int channel = static_cast<int>(m_decoder->get(cellID, "channel"));
        const int sca = static_cast<int>(m_decoder->get(cellID, "sca"));
        const bool fired = trigger ? ((*trigger)[hitIndex] > 0) : true;
        ++hitIndex;

        const double gainMpv = m_gainTable.at(slab, chip, channel);
        const double recoMpv = m_recoTable.at(slab, chip, channel);
        double pedestal = m_pedestals.at(slab, chip, channel, sca);
        if (std::isnan(pedestal)) pedestal = m_pedestalFallback.value();

        float adcHighPedsub = 0.0f;
        float adcLowPedsub = 0.0f;
        float energyMip = 0.0f;

        if (std::isnan(gainMpv) || std::isnan(recoMpv)) {
          // Not calibrated: EventBuilder::buildHit gives such a channel an energy
          // of 0 rather than a guess, and so does this.
          ++nUncalibrated;
        } else {
          // 1. the charge, in ADC, as the channel would really have given it.
          //    gainMpv is the channel's RELATIVE response (table / median), so
          //    the absolute scale is AdcPerMip and nothing else.
          const double idealPedsub = hit.getEnergy() * m_adcPerMip.value() * gainMpv;

          // 2. the readout: the high gain saturates, the low gain does not.
          //    The high gain does not clip at a single value, it rolls over:
          //    S(q) = q / (1 + (q/A)^n)^(1/n), linear well below A, saturating
          //    at A, with n setting how sharp the knee is (n -> inf is a hard
          //    clip).  A and n are the DYNAMIC RANGE of this threshold set and
          //    are fitted per set from the electron data's own hit_lg-vs-hit_hg
          //    relation (analysis/compare_gains.py), where the linear low gain
          //    gives the true charge and the high gain shows what it did with it.
          const double high = rollOver(idealPedsub);
          // The low gain is not a line in the data: around the anchor there is
          // a 3.4-3.9 LG-ADC pedestal width of its own, and a ~5 % dispersion
          // of the gain ratio from channel to channel (analysis/compare_gains.py
          // measures both, `lg_noise_adc` / `lg_gain_spread` in
          // digi_calibration.yml).  Both are symmetric, so they widen the
          // energy of the hits reconstructed from the low gain (above
          // SaturationAdc) without moving its mean; the spread is drawn per hit,
          // which gives a spectrum the same marginal as a per-channel factor.
          double low = m_gainRatio.value() * idealPedsub + m_gainIntercept.value();
          if (m_lowGainNoiseAdc.value() > 0.0 || m_lowGainSpread.value() > 0.0) {
            const auto seed = static_cast<std::uint64_t>(m_randomSeed.value())
                + 0x9e3779b97f4a7c15ULL * static_cast<std::uint64_t>(evtNum + 1)
                + static_cast<std::uint64_t>(hitIndex);
            std::mt19937_64 rng(seed);
            std::normal_distribution<double> gauss(0.0, 1.0);
            low = m_gainRatio.value() * idealPedsub * (1.0 + m_lowGainSpread.value() * gauss(rng))
                + m_gainIntercept.value() + m_lowGainNoiseAdc.value() * gauss(rng);
          }
          double highQ = high;
          if (m_quantise.value()) {
            highQ = std::round(high + pedestal) - pedestal;
            low = std::round(low + pedestal) - pedestal;
          }
          adcHighPedsub = static_cast<float>(highQ);
          adcLowPedsub = static_cast<float>(low);

          // 3. the reconstruction, port of EventBuilder::buildHit
          const bool saturated = highQ >= m_saturationAdc.value();
          if (saturated) {
            ++nSaturated;
            const double highEquivalent =
                (low - m_gainIntercept.value()) / m_gainRatio.value();
            energyMip = static_cast<float>(highEquivalent / recoMpv);
          } else {
            energyMip = static_cast<float>(highQ / recoMpv);
          }
        }

        auto nh = output->create();
        nh.setCellID(cellID);
        nh.setEnergy(energyMip);
        nh.setPosition(hit.getPosition());
        for (const auto& contrib : hit.getContributions()) {
          nh.addToContributions(contrib);
        }
        adcHighOut->push_back(adcHighPedsub);
        adcLowOut->push_back(adcLowPedsub);

        // The chip-level decision.  Without a trigger collection every hit is
        // kept (the cell-level selection already happened upstream).  With one:
        // a cell on a chip that fired is a hit if it fired itself or if its
        // sampled high gain is above the threshold; anything else is written
        // (the collections stay 1:1) but flagged out, and sim_to_ecal_tree
        // drops it.
        int kept = 1;
        if (trigger) {
          const bool onFiredChip = chipFired.count({slab, chip}) > 0;
          kept = (onFiredChip && (fired || adcHighPedsub > m_adcHitThreshold.value())) ? 1 : 0;
          if (kept && !fired) ++nRecovered;
          if (!kept) ++nDropped;
        }
        keptOut->push_back(kept);
      }

      if (doPrint) {
        debug() << "AdcDigitizer: " << input->size() << " hits, "
                << nUncalibrated << " uncalibrated (energy 0), "
                << nSaturated << " over the saturation threshold"
                << (trigger ? (", chip selection: " + std::to_string(nRecovered) +
                               " kept without their own trigger, " + std::to_string(nDropped) +
                               " dropped")
                            : std::string())
                << endmsg;
      }
      return StatusCode::SUCCESS;
    } catch (const std::exception& e) {
      error() << "[AdcDigitizer] Exception in execute(): " << e.what() << endmsg;
      return StatusCode::FAILURE;
    }
  }

  StatusCode finalize() override {
    m_inputHandle.reset();
    m_outputHandle.reset();
    m_adcHighHandle.reset();
    m_adcLowHandle.reset();
    m_keptHandle.reset();
    m_triggerHandle.reset();
    m_decoder.reset();
    return Gaudi::Algorithm::finalize();
  }

private:
  bool loadTables() {
    const std::string gainPath = tablePath("mips", m_gainShapeThreshold.value(), "MIP_");
    const std::string recoPath = tablePath("mips", m_calibThreshold.value(), "MIP_");
    const std::string pedPath = tablePath("pedestals", m_calibThreshold.value(), "Pedestal_");
    std::string err;

    if (gainPath.empty() || !siwecal::readMipTable(gainPath, m_gainTable, err)) {
      error() << "[AdcDigitizer] gain-shape table: "
              << (err.empty() ? "not found" : err) << endmsg;
      return false;
    }
    if (!normaliseGainShape()) return false;
    if (recoPath.empty() || !siwecal::readMipTable(recoPath, m_recoTable, err)) {
      error() << "[AdcDigitizer] reconstruction table: "
              << (err.empty() ? "not found" : err) << endmsg;
      return false;
    }
    if (pedPath.empty() || !siwecal::readPedestalTable(pedPath, m_pedestals, err)) {
      // The pedestal has exactly one job here -- it sets the grid the ADC is
      // rounded onto -- because everything downstream works pedestal-subtracted.
      // With an integer fallback the grid is the same as with the real table, so
      // this costs nothing measurable; it is an info, not a missing input.
      info() << "[AdcDigitizer] no pedestal table under " << m_calibDir.value()
             << " (" << (err.empty() ? "not found" : err) << "); using "
             << m_pedestalFallback.value()
             << " ADC, which only sets the rounding grid" << endmsg;
    }
    if (m_gainTable.nCalibrated == 0 || m_recoTable.nCalibrated == 0) {
      error() << "[AdcDigitizer] a calibration table came back empty." << endmsg;
      return false;
    }
    return true;
  }

  /// The high-gain roll-over: identity for q << A, saturating at A.
  double rollOver(double q) const {
    const double a = m_adcHighMax.value();
    const double n = m_saturationOrder.value();
    if (q <= 0.0 || a <= 0.0) return q;
    if (n <= 0.0) return std::min(q, a);            // hard clip, the old model
    return q / std::pow(1.0 + std::pow(q / a, n), 1.0 / n);
  }

  /// Divide the gain table by its own median, so it carries only the relative
  /// channel-to-channel variation and AdcPerMip carries the scale.
  static double medianOf(std::vector<double> values) {
    const std::size_t mid = values.size() / 2;
    std::nth_element(values.begin(), values.begin() + mid, values.end());
    return values[mid];
  }

  // The gain table is a SHAPE: channel-to-channel variation only, the absolute
  // scale is AdcPerMip.  "global" divides every channel by the detector median.
  // "per-slab" divides each slab by its own median, so a slab whose table MPV is
  // off for a reason that is not gain never enters as gain: slab 12 of the 2026
  // stack ran at a higher discriminator (DAC 243) and its truncated MIP spectra
  // put it at 1.2-1.9x the others depending on which table is asked -- a real
  // gain would not depend on the threshold set -- and layer 14's 650 um sensor
  // already collects 1.34x the charge in the simulation itself, so the table
  // must not add it again.  Per slab, the data reconstruction does the same
  // thing implicitly: it divides every channel by its own table.
  bool normaliseGainShape() {
    std::vector<double> all;
    all.reserve(m_gainTable.mpv.size());
    for (double v : m_gainTable.mpv) {
      if (std::isfinite(v)) all.push_back(v);
    }
    if (all.empty()) {
      error() << "[AdcDigitizer] gain-shape table has no calibrated channel."
              << endmsg;
      return false;
    }
    const double global = medianOf(all);
    if (!(global > 0.0)) {
      error() << "[AdcDigitizer] gain-shape median is " << global << endmsg;
      return false;
    }
    const std::string mode = m_gainShapeNormalisation.value();
    if (mode == "global") {
      for (double& v : m_gainTable.mpv) v /= global;
      info() << "[AdcDigitizer] gain shape normalised by its median "
             << global << " ADC/MIP" << endmsg;
      return true;
    }
    if (mode != "per-slab") {
      error() << "[AdcDigitizer] GainShapeNormalisation must be 'global' or 'per-slab', got '"
              << mode << "'" << endmsg;
      return false;
    }
    info() << "[AdcDigitizer] gain shape normalised per slab (detector median "
           << global << " ADC/MIP); slab medians / global:";
    for (int slab = 0; slab < siwecal::kSlabs; ++slab) {
      std::vector<double> vals;
      const std::size_t first = static_cast<std::size_t>(slab) * siwecal::kChips *
                                siwecal::kChannels;
      const std::size_t last = first + siwecal::kChips * siwecal::kChannels;
      for (std::size_t i = first; i < last; ++i) {
        if (std::isfinite(m_gainTable.mpv[i])) vals.push_back(m_gainTable.mpv[i]);
      }
      const double median = vals.empty() ? global : medianOf(vals);
      for (std::size_t i = first; i < last; ++i) m_gainTable.mpv[i] /= median;
      info() << " " << slab << ":" << std::fixed << std::setprecision(2) << median / global;
    }
    info() << endmsg;
    return true;
  }

  /// <CalibDir>/<kind>/<threshold>/<prefix>*_<gain>.txt, resolved by globbing the
  /// directory the same way ChannelMapper does.
  std::string tablePath(const std::string& kind, const std::string& threshold,
                        const std::string& prefix) const {
    namespace fs = std::filesystem;
    const fs::path dir = fs::path(m_calibDir.value()) / kind / threshold;
    std::error_code ec;
    if (!fs::is_directory(dir, ec)) return {};
    for (const auto& entry : fs::directory_iterator(dir, ec)) {
      const std::string name = entry.path().filename().string();
      if (name.rfind(prefix, 0) == 0 &&
          name.find(m_calibGain.value()) != std::string::npos) {
        return entry.path().string();
      }
    }
    return {};
  }

  Gaudi::Property<std::string> m_inputName{
      this, "InputCollection", "SiPadHitsRealMapped",
      "Input hits, already in test-beam CellID format, energy in MIP"};
  Gaudi::Property<std::string> m_outputName{
      this, "OutputCollection", "SiPadHitsRealAdc",
      "Output hits whose energy is the RECONSTRUCTED MIP, i.e. what the same "
      "calibration the data uses makes of the simulated ADC"};
  Gaudi::Property<std::string> m_adcHighName{
      this, "AdcHighCollection", "SiPadHitsRealAdcHigh",
      "Pedestal-subtracted high-gain ADC per hit (the data's hit_hg)"};
  Gaudi::Property<std::string> m_adcLowName{
      this, "AdcLowCollection", "SiPadHitsRealAdcLow",
      "Pedestal-subtracted low-gain ADC per hit (the data's hit_lg)"};
  Gaudi::Property<std::string> m_keptName{
      this, "KeptCollection", "SiPadHitsRealAdcKept",
      "Output UserDataCollection<int>, 1:1 with the output hits: 1 if the hit "
      "survives the chip-level selection (always 1 without TriggerCollection)"};
  Gaudi::Property<std::string> m_triggerName{
      this, "TriggerCollection", "",
      "Input UserDataCollection<int> of per-cell trigger flags (RealDigitizer's "
      "DigitizedTriggerCollection, carried 1:1 through the flipper and the mapper). "
      "Empty: no chip-level selection, every hit is kept"};
  Gaudi::Property<double> m_adcHitThreshold{
      this, "AdcHitThreshold", 30.0,
      "Chip-level selection: pedestal-subtracted high gain above which a cell that did "
      "not fire is kept on a chip that did (EcalEventBuilder's AdcHitThreshold)"};
  Gaudi::Property<std::string> m_calibDir{
      this, "CalibDir", "masking_info/calibration/MuonCalib_gaudi",
      "Root of the MuonCalib_gaudi tree"};
  Gaudi::Property<std::string> m_calibThreshold{
      this, "CalibThreshold", "th230",
      "Threshold set whose MIP table RECONSTRUCTS the hit -- the one the data "
      "run was reconstructed with"};
  Gaudi::Property<double> m_adcPerMip{
      this, "AdcPerMip", 18.564,
      "ADC per MIP: the absolute charge scale, fitted per threshold set against "
      "test-beam data at a reference energy by analysis/fit_adc_scale.py and "
      "stored in mappings/digi_calibration.yml. The default is only a starting "
      "point for the first iteration of that fit"};
  Gaudi::Property<std::string> m_gainShapeThreshold{
      this, "GainShapeThreshold", "th210",
      "Threshold set whose MIP table supplies the channel-to-channel VARIATION "
      "of the gain (normalised to its own median). Only the shape is taken from "
      "it, never the scale -- that is AdcPerMip"};
  Gaudi::Property<std::string> m_gainShapeNormalisation{
      this, "GainShapeNormalisation", "global",
      "'global': the gain-shape table is divided by the detector median; 'per-slab': each "
      "slab by its own median, so a slab-level table offset (a truncated MIP fit, a thicker "
      "sensor already in the geometry) does not enter as gain"};
  Gaudi::Property<std::string> m_calibGain{
      this, "CalibGain", "highgain", "highgain or lowgain"};
  Gaudi::Property<double> m_adcHighMax{
      this, "AdcHighMax", 2100.0,
      "A of the roll-over S(q) = q/(1+(q/A)^n)^(1/n): the ceiling the high gain "
      "saturates at [pedestal-subtracted ADC]. Fitted per threshold set from the "
      "electron data's hit_lg-vs-hit_hg relation (analysis/compare_gains.py)"};
  Gaudi::Property<double> m_saturationOrder{
      this, "SaturationOrder", 8.0,
      "n of the same roll-over: how sharp the knee is. 0 means a hard clip at "
      "AdcHighMax. Fitted together with AdcHighMax"};
  Gaudi::Property<double> m_saturationAdc{
      this, "SaturationAdc", 1500.0,
      "Pedestal-subtracted high-gain ADC at which the reconstruction switches to "
      "the low gain (EcalEventBuilder.AdcSaturationThreshold)"};
  Gaudi::Property<double> m_gainRatio{
      this, "GainRatio", 0.0961, "k of adc_low - ped_lg = k*(adc_high - ped_hg) + c"};
  Gaudi::Property<double> m_gainIntercept{
      this, "GainIntercept", 1.60, "c of the same line"};
  Gaudi::Property<double> m_lowGainNoiseAdc{
      this, "LowGainNoiseAdc", 0.0,
      "Gaussian noise added to the low gain [LG ADC], the width of the data's "
      "hit_lg about the anchor line for small hits (digi_calibration.yml "
      "lg_noise_adc). 0 keeps the line"};
  Gaudi::Property<double> m_lowGainSpread{
      this, "LowGainSpread", 0.0,
      "Relative Gaussian spread of the low-gain slope per hit (digi_calibration.yml "
      "lg_gain_spread): the part of the data's scatter that grows with the "
      "amplitude, i.e. the channel-to-channel dispersion of the gain ratio"};
  Gaudi::Property<unsigned long long> m_randomSeed{
      this, "RandomSeed", 7919ULL, "Base seed; each hit draws from seed+event+index"};
  Gaudi::Property<double> m_pedestalFallback{
      this, "PedestalFallback", 245.0,
      "Pedestal used where the table has none [ADC]"};
  Gaudi::Property<bool> m_quantise{
      this, "Quantise", true, "Round the ADC to integers, as the chip does"};
  Gaudi::Property<std::string> m_bitField{
      this, "BitField", "system:8,slab:8,chip:16,channel:8,sca:8",
      "Test-beam CellID bitfield"};
  Gaudi::Property<int> m_debugFreq{this, "DebugFrequency", 500, ""};

  mutable std::unique_ptr<k4FWCore::DataHandle<edm4hep::SimCalorimeterHitCollection>> m_inputHandle;
  mutable std::unique_ptr<k4FWCore::DataHandle<edm4hep::SimCalorimeterHitCollection>> m_outputHandle;
  mutable std::unique_ptr<k4FWCore::DataHandle<podio::UserDataCollection<float>>> m_adcHighHandle;
  mutable std::unique_ptr<k4FWCore::DataHandle<podio::UserDataCollection<float>>> m_adcLowHandle;
  mutable std::unique_ptr<k4FWCore::DataHandle<podio::UserDataCollection<int>>> m_keptHandle;
  mutable std::unique_ptr<k4FWCore::DataHandle<podio::UserDataCollection<int>>> m_triggerHandle;
  std::unique_ptr<dd4hep::DDSegmentation::BitFieldCoder> m_decoder;
  siwecal::MipTable m_gainTable;
  siwecal::MipTable m_recoTable;
  siwecal::PedestalTable m_pedestals;
  mutable std::atomic<long long> m_eventCount{0};
};

DECLARE_COMPONENT(AdcDigitizer)

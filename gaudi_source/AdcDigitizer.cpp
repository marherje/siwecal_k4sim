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
#include <fstream>
#include <sstream>
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
//   1. MIP -> ADC with `AdcPerMip`, the ADC a channel gives per MIP: the only
//      absolute scale in the chain, and ONE number for every threshold set,
//      because the preamplifier gain does not depend on the trigger DAC.  It is
//      measured on muons where the discriminator sits below the MIP peak (th210:
//      spectrum shape, analysis/muon_gain_scan.py, and the tag-and-probe layer
//      efficiency of muon tracks), 19.5 +- 0.5 ADC/MIP (2026-09-22).  A per-set
//      muon anchor was tried and is not a measurement above the MIP: the
//      surviving Landau tail is scale-free, and the scan's chi2 without the
//      simulation's own variance is its noise floor.  Neither is a shower
//      anchor: it absorbs whatever the simulation gets wrong about the shower.
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
// 4. (OccupancySlope != 0) the chip's baseline sag under its own load: the
//    other hits of the same chip pull the common baseline by B = k * Q, Q = the
//    HG-ADC charge of the chip's other hits.  Measured on pedestals
//    (k = -0.95e-3) and fitted on the electron chip multiplicity
//    (k = -1.5e-3, sim_tuning 2.1b).  A cell read in the high gain loses |B| in
//    amplitude and energy, and its trigger survives with P(a + B) / P(a), the
//    per-slab erf turn-on at its pre-sag amplitude a.  The low gain moves by
//    GainRatio * B, which is negligible, and is left alone.  Q is taken before
//    the sag (the shift is small against the charge it comes from), from the
//    hits the chip would keep, at or above OccupancyMinMip.
//    With OccupancySagTable (2026-10-29, replaces k * Q): B is read from a map
//    MEASURED on the electron runs of the threshold set, B(Q_chip, Q_slab), with
//    Q_chip as above and Q_slab = the HG-ADC charge of the other chips of the
//    same slab: the baseline also drops with the load of the whole slab, and
//    saturates (analysis/cm_sag_collect.py, sag_table.py).  Bilinear in
//    (ln(1+Q_chip), ln(1+Q_slab)) between the map's bin centres, clamped at its
//    edges.
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

      if (!m_occupancySagTable.value().empty()) {
        if (!loadSagTable(m_occupancySagTable.value())) return StatusCode::FAILURE;
        if (m_occupancySlope.value() != 0.0) {
          error() << "[AdcDigitizer] OccupancySlope and OccupancySagTable are exclusive" << endmsg;
          return StatusCode::FAILURE;
        }
      }
      if (sagOn()) {
        if (m_occTriggerAdc.value().empty() || m_occTriggerSigmaAdc.value().empty()) {
          error() << "[AdcDigitizer] OccupancySlope needs OccupancyTriggerAdc and OccupancyTriggerSigmaAdc "
                     "(the per-slab turn-on in HG ADC)" << endmsg;
          return StatusCode::FAILURE;
        }
        info() << "[AdcDigitizer] occupancy baseline sag ON: "
               << (m_sagQc.empty() ? "B = " + std::to_string(m_occupancySlope.value()) + " x Q"
                                   : "B(Q_chip, Q_slab) map " + m_occupancySagTable.value() + " (" +
                                         std::to_string(m_sagQc.size()) + " x " + std::to_string(m_sagQs.size()) +
                                         ", B " + std::to_string(m_sagMin) + " .. " + std::to_string(m_sagMax) + ")")
               << " [HG ADC], Q = the chip's other kept hits >= " << m_occupancyMinMip.value()
               << " MIP; trigger survival P(a+B)/P(a) with " << m_occTriggerAdc.value().size()
               << " per-slab turn-ons" << endmsg;
      }

      m_decoder = std::make_unique<dd4hep::DDSegmentation::BitFieldCoder>(m_bitField.value());

      info() << "[AdcDigitizer] AdcPerMip=" << m_adcPerMip.value()
             << " ADC/MIP (shape from " << m_gainShapeThreshold.value() << ", "
             << m_gainTable.nCalibrated << " channels), reconstruction table "
             << (m_recoTableThreshold.value().empty() ? m_calibThreshold.value()
                                                      : m_recoTableThreshold.value())
             << " (" << m_recoTable.nCalibrated
             << " channels), pedestals " << m_pedestals.nCalibrated << " entries"
             << (m_pedestalNoise.value() ? ", per-cell pedestal noise ON" : "") << endmsg;
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

      // Pass 1: the readout of every cell.  The writing waits for pass 2, because the
      // chip's baseline sag (step 4) needs the whole chip first.
      const std::size_t nIn = input->size();
      std::vector<int> vSlab(nIn), vChip(nIn), vFired(nIn), vCalib(nIn), vSat(nIn);
      std::vector<double> vHigh(nIn), vRecoMpv(nIn);
      std::vector<float> vHighOut(nIn), vLowOut(nIn), vEnergy(nIn);

      std::size_t hitIndex = 0;
      for (const auto& hit : *input) {
        const std::size_t ih = hitIndex;
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
          double high = rollOver(idealPedsub);
          // Electronic noise of the cell, from the pedestal table's own fitted width
          // (per channel and SCA), after the roll-over: it is the readout's, not the
          // charge's.  The pedestal widths of the fixed-SCA-pairing tables (median
          // 1.37 ADC at th210) are the measured noise; a scalar SlowNoiseMIP in
          // RealDigitizer stands in for it when this is off.
          if (m_pedestalNoise.value()) {
            // Per (channel, SCA) when the SCA is real, else the channel's mean over its SCAs.
            double w = m_pedestalNoisePerSca.value() ? m_pedestals.widthAt(slab, chip, channel, sca)
                                                     : m_pedestals.meanWidth(slab, chip, channel);
            if (std::isnan(w)) w = m_pedestalNoiseFallback.value();
            if (w > 0.0) {
              const auto seed = static_cast<std::uint64_t>(m_randomSeed.value()) * 31ULL
                  + 0xbf58476d1ce4e5b9ULL * static_cast<std::uint64_t>(evtNum + 1)
                  + static_cast<std::uint64_t>(hitIndex);
              std::mt19937_64 rng(seed);
              std::normal_distribution<double> gauss(0.0, 1.0);
              high += w * gauss(rng);
            }
          }
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
          vHigh[ih] = highQ;
          vRecoMpv[ih] = recoMpv;
          vCalib[ih] = 1;

          // 3. the reconstruction, port of EventBuilder::buildHit
          const bool saturated = highQ >= m_saturationAdc.value();
          vSat[ih] = saturated ? 1 : 0;
          if (saturated) {
            ++nSaturated;
            const double highEquivalent =
                (low - m_gainIntercept.value()) / m_gainRatio.value();
            energyMip = static_cast<float>(highEquivalent / recoMpv);
          } else {
            energyMip = static_cast<float>(highQ / recoMpv);
          }
        }

        vSlab[ih] = slab;
        vChip[ih] = chip;
        vFired[ih] = fired ? 1 : 0;
        vHighOut[ih] = adcHighPedsub;
        vLowOut[ih] = adcLowPedsub;
        vEnergy[ih] = energyMip;
      }

      // The chip-level decision.  Without a trigger collection every hit is
      // kept (the cell-level selection already happened upstream).  With one:
      // a cell on a chip that fired is a hit if it fired itself or if its
      // sampled high gain is above the threshold; anything else is written
      // (the collections stay 1:1) but flagged out, and sim_to_ecal_tree
      // drops it.
      auto keepRule = [&](std::size_t i) {
        if (!trigger) return vFired[i] > 0;
        const bool onFiredChip = chipFired.count({vSlab[i], vChip[i]}) > 0;
        return onFiredChip && (vFired[i] > 0 || vHighOut[i] > m_adcHitThreshold.value());
      };

      // 4. the chip's baseline sag, B = k * Q or the measured B(Q_chip, Q_slab) (see the header).
      int nOccLost = 0;
      double occShiftSum = 0.0;
      if (sagOn()) {
        std::map<std::pair<int, int>, double> chipQ;
        std::map<int, double> slabQ;
        std::vector<double> ownQ(nIn, 0.0);
        for (std::size_t i = 0; i < nIn; ++i) {
          if (!keepRule(i) || vEnergy[i] < m_occupancyMinMip.value()) continue;
          ownQ[i] = vEnergy[i] * m_adcPerMip.value();
          chipQ[{vSlab[i], vChip[i]}] += ownQ[i];
          slabQ[vSlab[i]] += ownQ[i];
        }
        for (std::size_t i = 0; i < nIn; ++i) {
          const auto it = chipQ.find({vSlab[i], vChip[i]});
          const double qChipAll = it == chipQ.end() ? 0.0 : it->second;
          const double q = qChipAll - ownQ[i];
          double b;
          if (m_sagQc.empty()) {
            b = m_occupancySlope.value() * std::max(q, 0.0);
          } else {
            const auto is = slabQ.find(vSlab[i]);
            const double qSlab = (is == slabQ.end() ? 0.0 : is->second) - qChipAll;   // the OTHER chips
            b = sagFromTable(std::max(q, 0.0), std::max(qSlab, 0.0));
          }
          if (b == 0.0 || !vCalib[i]) continue;
          // the trigger, at the amplitude before the sag
          if (vFired[i]) {
            const double mu = perSlab(m_occTriggerAdc.value(), vSlab[i]);
            const double sg = perSlab(m_occTriggerSigmaAdc.value(), vSlab[i]);
            const double a = vHigh[i];
            const double p0 = turnOn(a, mu, sg), p1 = turnOn(a + b, mu, sg);
            const double survive = p0 > 1e-9 ? std::clamp(p1 / p0, 0.0, 1.0) : 1.0;
            const auto seed = static_cast<std::uint64_t>(m_randomSeed.value()) * 131ULL
                + 0x94d049bb133111ebULL * static_cast<std::uint64_t>(evtNum + 1)
                + static_cast<std::uint64_t>(i);
            std::mt19937_64 rng(seed);
            if (std::uniform_real_distribution<double>(0.0, 1.0)(rng) >= survive) {
              vFired[i] = 0;
              ++nOccLost;
            }
          }
          // the amplitude and energy of a cell read in the high gain
          if (!vSat[i]) {
            vHigh[i] += b;
            vHighOut[i] = static_cast<float>(vHigh[i]);
            vEnergy[i] = static_cast<float>(vHigh[i] / vRecoMpv[i]);
          }
          occShiftSum += b;
        }
        if (trigger) {
          chipFired.clear();
          for (std::size_t i = 0; i < nIn; ++i)
            if (vFired[i]) chipFired[{vSlab[i], vChip[i]}] = true;
        }
      }

      // Pass 2: write, 1:1 with the input.
      hitIndex = 0;
      for (const auto& hit : *input) {
        const std::size_t i = hitIndex++;
        auto nh = output->create();
        nh.setCellID(hit.getCellID());
        nh.setEnergy(vEnergy[i]);
        nh.setPosition(hit.getPosition());
        for (const auto& contrib : hit.getContributions()) {
          nh.addToContributions(contrib);
        }
        adcHighOut->push_back(vHighOut[i]);
        adcLowOut->push_back(vLowOut[i]);
        const int kept = keepRule(i) ? 1 : 0;
        if (trigger) {
          if (kept && !vFired[i]) ++nRecovered;
          if (!kept) ++nDropped;
        }
        keptOut->push_back(kept);
      }
      m_nOccLost += nOccLost;
      m_nOccEvents += 1;

      if (doPrint) {
        debug() << "AdcDigitizer: " << input->size() << " hits, "
                << nUncalibrated << " uncalibrated (energy 0), "
                << nSaturated << " over the saturation threshold"
                << (trigger ? (", chip selection: " + std::to_string(nRecovered) +
                               " kept without their own trigger, " + std::to_string(nDropped) +
                               " dropped")
                            : std::string())
                << (sagOn()
                        ? ", occupancy: " + std::to_string(nOccLost) + " triggers lost, summed shift " +
                              std::to_string(occShiftSum) + " ADC"
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
    if (sagOn() && m_nOccEvents.load() > 0) {
      info() << "[AdcDigitizer] occupancy baseline sag: " << m_nOccLost.load() << " triggers lost in "
             << m_nOccEvents.load() << " events ("
             << static_cast<double>(m_nOccLost.load()) / m_nOccEvents.load() << " per event)" << endmsg;
    }
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
    // The table that turns the ADC back into MIP: the run's own set by default,
    // or RecoTableThreshold when the data were reconstructed with another set's
    // table (the Reconstruction_adc_th210 campaign converts every run with the
    // th210 table, the only one fitted with the discriminator below the MIP
    // peak).  Pedestals and, in ChannelMapper, the mask stay with CalibThreshold.
    const std::string recoTh = m_recoTableThreshold.value().empty() ? m_calibThreshold.value()
                                                                    : m_recoTableThreshold.value();
    const std::string recoPath = tablePath("mips", recoTh, "MIP_");
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

  /// Erf turn-on of the trigger at HG amplitude a [ADC].
  bool sagOn() const { return m_occupancySlope.value() != 0.0 || !m_sagQc.empty(); }

  // OccupancySagTable: '#' comments, 'qchip <centres>', 'qslab <centres>', one row of B per Q_chip centre.
  bool loadSagTable(const std::string& path) {
    std::ifstream in(path);
    if (!in) {
      error() << "[AdcDigitizer] cannot open OccupancySagTable " << path << endmsg;
      return false;
    }
    std::string line;
    std::vector<double> qc, qs, b;
    while (std::getline(in, line)) {
      if (line.empty() || line[0] == '#') continue;
      std::istringstream ss(line);
      std::string tag;
      ss >> tag;
      std::vector<double>* dst = tag == "qchip" ? &qc : tag == "qslab" ? &qs : nullptr;
      if (!dst) {
        std::istringstream row(line);
        double v;
        while (row >> v) b.push_back(v);
        continue;
      }
      double v;
      while (ss >> v) dst->push_back(std::log1p(v));
    }
    if (qc.size() < 2 || qs.size() < 2 || b.size() != qc.size() * qs.size()) {
      error() << "[AdcDigitizer] OccupancySagTable " << path << ": " << qc.size() << " x " << qs.size()
              << " centres but " << b.size() << " values" << endmsg;
      return false;
    }
    m_sagQc = qc;
    m_sagQs = qs;
    m_sagB = b;
    m_sagMin = *std::min_element(b.begin(), b.end());
    m_sagMax = *std::max_element(b.begin(), b.end());
    return true;
  }

  // Bilinear in (ln(1+Q_chip), ln(1+Q_slab)), clamped at the edges (as analysis/sag_table.py interp).
  double sagFromTable(double qChip, double qSlab) const {
    auto locate = [](const std::vector<double>& x, double v, std::size_t& i, double& t) {
      v = std::clamp(v, x.front(), x.back());
      i = std::min<std::size_t>(std::upper_bound(x.begin(), x.end(), v) - x.begin(), x.size() - 1);
      i = i == 0 ? 0 : i - 1;
      t = (v - x[i]) / (x[i + 1] - x[i]);
    };
    std::size_t i, j;
    double tu, tv;
    locate(m_sagQc, std::log1p(qChip), i, tu);
    locate(m_sagQs, std::log1p(qSlab), j, tv);
    const std::size_t n = m_sagQs.size();
    auto B = [&](std::size_t r, std::size_t c) { return m_sagB[r * n + c]; };
    return (1 - tu) * (1 - tv) * B(i, j) + tu * (1 - tv) * B(i + 1, j) + (1 - tu) * tv * B(i, j + 1) +
           tu * tv * B(i + 1, j + 1);
  }

  static double turnOn(double a, double mu, double sg) {
    return 0.5 * (1.0 + std::erf((a - mu) / (std::sqrt(2.0) * std::max(sg, 1e-6))));
  }

  /// v[slab], or its last entry for a slab past the end (a single value = every slab).
  static double perSlab(const std::vector<double>& v, int slab) {
    return v[std::min<std::size_t>(static_cast<std::size_t>(std::max(slab, 0)), v.size() - 1)];
  }

  /// The high-gain roll-over: identity for q << A, saturating at A.
  double rollOver(double q) const {
    const double a = m_adcHighMax.value();
    const double n = m_saturationOrder.value();
    if (q <= 0.0 || a <= 0.0) return q;
    if (n <= 0.0) return std::min(q, a);            // hard clip, the old model
    return q / std::pow(1.0 + std::pow(q / a, n), 1.0 / n);
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
    std::string err, log;
    if (!siwecal::normaliseGainShape(m_gainTable, m_gainShapeNormalisation.value(), m_gainShapeShrink.value(),
                                     m_gainShapeTrueSpread.value(), err, log)) {
      error() << "[AdcDigitizer] " << err << endmsg;
      return false;
    }
    std::istringstream lines(log);
    for (std::string line; std::getline(lines, line);) info() << "[AdcDigitizer] " << line << endmsg;
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
  Gaudi::Property<std::string> m_recoTableThreshold{
      this, "RecoTableThreshold", "",
      "Threshold set whose MIP table turns the ADC back into MIP, when it is not "
      "the run's own (empty = CalibThreshold). Must be the table the data were "
      "reconstructed with: th210 for the Reconstruction_adc_th210 campaign. The "
      "pedestals stay with CalibThreshold"};
  Gaudi::Property<double> m_adcPerMip{
      this, "AdcPerMip", 19.5,
      "ADC per MIP: the absolute charge scale, one value for every threshold set, "
      "measured on muons where the discriminator sits below the MIP peak "
      "(analysis/muon_gain_scan.py) and stored in mappings/digi_calibration.yml"};
  Gaudi::Property<std::string> m_gainShapeThreshold{
      this, "GainShapeThreshold", "th210",
      "Threshold set whose MIP table supplies the channel-to-channel VARIATION "
      "of the gain (normalised to its own median). Only the shape is taken from "
      "it, never the scale -- that is AdcPerMip"};
  Gaudi::Property<bool> m_gainShapeShrink{
      this, "GainShapeShrink", false,
      "Remove the MPV fit error from the gain shape's channel-to-channel spread "
      "(per-slab linear shrinkage toward the slab median, see siwecal::normaliseGainShape)"};
  Gaudi::Property<double> m_gainShapeTrueSpread{
      this, "GainShapeTrueSpread", 0.0,
      "With GainShapeShrink: the true relative channel-to-channel gain spread per slab. > 0 scales each "
      "slab's deviations to exactly this spread (measured as sqrt(cov) of two independent MIP tables, "
      "th210 x th220 fixed: 0.032); 0 estimates it from the table's own MPV errors"};
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
  Gaudi::Property<bool> m_pedestalNoise{
      this, "PedestalNoise", false,
      "Add each cell's electronic noise to the high gain: a Gaussian of the pedestal "
      "table's fitted width at (slab, chip, channel, sca). Use with RealDigitizer "
      "SlowNoiseMIP = 0, which otherwise adds a scalar noise before the ADC"};
  Gaudi::Property<bool> m_pedestalNoisePerSca{
      this, "PedestalNoisePerSca", false,
      "Take the noise of the hit's own SCA; false (default, the simulation's sca is always 0) "
      "averages the channel's widths over its fitted SCAs"};
  Gaudi::Property<double> m_pedestalNoiseFallback{
      this, "PedestalNoiseFallback", 1.37,
      "Noise [ADC] where the pedestal table has no width (median th210 width, fixed "
      "SCA pairing)"};
  Gaudi::Property<double> m_pedestalFallback{
      this, "PedestalFallback", 245.0,
      "Pedestal used where the table has none [ADC]"};
  Gaudi::Property<bool> m_quantise{
      this, "Quantise", true, "Round the ADC to integers, as the chip does"};
  Gaudi::Property<std::string> m_bitField{
      this, "BitField", "system:8,slab:8,chip:16,channel:8,sca:8",
      "Test-beam CellID bitfield"};
  Gaudi::Property<double> m_occupancySlope{
      this, "OccupancySlope", 0.0,
      "k of the chip baseline sag B = k * Q [HG ADC per HG ADC], Q = HG-ADC charge of the chip's other hits. "
      "0 = off. Fitted on the electron chip multiplicity: -1.5e-3 (digi_calibration.yml occupancy_slope)"};
  Gaudi::Property<std::string> m_occupancySagTable{
      this, "OccupancySagTable", "",
      "Measured chip baseline sag map B(Q_chip, Q_slab) [HG ADC] (digi_calibration.yml occupancy_sag_table); "
      "replaces OccupancySlope, which must then be 0. Empty = off"};
  std::vector<double> m_sagQc, m_sagQs, m_sagB;   // ln(1+Q) centres and the row-major map
  double m_sagMin = 0.0, m_sagMax = 0.0;
  Gaudi::Property<double> m_occupancyMinMip{
      this, "OccupancyMinMip", 0.5,
      "Hits below this energy [MIP] do not count in Q (the 0.5 MIP hit cut the fit was made with)"};
  Gaudi::Property<std::vector<double>> m_occTriggerAdc{
      this, "OccupancyTriggerAdc", {},
      "Per-slab trigger turn-on mean [HG ADC] for the survival P(a+B)/P(a) (threshold_adc / slab_overrides)"};
  Gaudi::Property<std::vector<double>> m_occTriggerSigmaAdc{
      this, "OccupancyTriggerSigmaAdc", {},
      "Per-slab trigger turn-on width [HG ADC] (sigma_mip * adc_per_mip / slab_overrides sigma_adc)"};
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
  mutable std::atomic<long long> m_nOccLost{0};
  mutable std::atomic<long long> m_nOccEvents{0};
};

DECLARE_COMPONENT(AdcDigitizer)

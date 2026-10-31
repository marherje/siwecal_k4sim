#include "Gaudi/Algorithm.h"
#include "GaudiKernel/MsgStream.h"
#include "k4FWCore/DataHandle.h"
#include "CellShaping.hh"
#include "CalibTables.hh"
#include "PadMap.hh"
#include "edm4hep/MCParticle.h"
#include "edm4hep/SimCalorimeterHitCollection.h"
#include "podio/UserDataCollection.h"
#include "DDSegmentation/BitFieldCoder.h"

#include <atomic>
#include <cmath>
#include <fstream>
#include <map>
#include <sstream>
#include <tuple>
#include <cstdint>
#include <memory>
#include <numeric>
#include <random>
#include <vector>

// Realistic SiPad digitizer.
//
// DigitizationMode:
//   "simple" — identical to BasicDigitizer: apply an energy threshold cut.
//   "real"   — contribution-level CR-RC shaping (include/CellShaping.hh).  Each
//              cell's CaloHitContributions become (energy, time) steps; the fast
//              channel gives the trigger time at its first threshold crossing and
//              the slow channel is sampled at triggerTime + DelayNs.  No full
//              waveform is built (that is the wave_scan path, for plotting only).
//
// InputEnergyUnit:
//   "MIP" — input hit energy is already in MIP units.
//   "GeV" — input hit energy is converted with the MIP calibration below.
//
// MIP calibration, mirroring GeV2MIPConversion:
//   MIPValues (per-layer vector) wins over MIPValue (scalar) when non-empty.
//   Per-layer mode needs BitField to decode 'layer' out of the CellID; NLayers
//   validates the vector length (0 skips the check).
//
// Mode "real" must run BEFORE any GeV->MIP conversion: GeV2MIPConversion (like
// BasicDigitizer, DetectorFlipper and ChannelMapper) creates fresh hits without
// copying Contributions, and the shaping needs them.  It normalises to MIP
// itself, so GeV2MIPConversion is redundant on this branch.
//
// HitEnergyContent decides what lands in the output hit's energy, which is what
// DetectorFlipper, ChannelMapper and analysis/sim_to_ecal_tree.py all read:
//   "digitized" — the shaped slow-sample amplitude [MIP].  This is what makes
//                 the digitised energy reach the analysis.
//   "input"     — the original input hit energy, untouched.
// Either way the shaped amplitude and the trigger time are also written to the
// parallel UserDataCollections named by DigitizedEnergyCollection /
// DigitizedTimeCollection, index-aligned with the output hits.
class RealDigitizer : public Gaudi::Algorithm {
public:
  RealDigitizer(const std::string& name, ISvcLocator* svcLoc)
      : Gaudi::Algorithm(name, svcLoc) {}

  StatusCode initialize() override {
    try {
      m_inputHandle  = std::make_unique<k4FWCore::DataHandle<edm4hep::SimCalorimeterHitCollection>>(
          m_inputName.value(),  Gaudi::DataHandle::Reader, this);
      m_outputHandle = std::make_unique<k4FWCore::DataHandle<edm4hep::SimCalorimeterHitCollection>>(
          m_outputName.value(), Gaudi::DataHandle::Writer, this);

      const std::string& mode = m_mode.value();
      if (mode != "simple" && mode != "real") {
        error() << "[RealDigitizer] Unknown DigitizationMode '" << mode
                << "'. Allowed values: 'simple', 'real'." << endmsg;
        return StatusCode::FAILURE;
      }
      const std::string& unit = m_inputEnergyUnit.value();
      if (unit != "MIP" && unit != "GeV") {
        error() << "[RealDigitizer] Unknown InputEnergyUnit '" << unit
                << "'. Allowed values: 'MIP', 'GeV'." << endmsg;
        return StatusCode::FAILURE;
      }
      const std::string& hitEnergy = m_hitEnergyContent.value();
      if (hitEnergy != "digitized" && hitEnergy != "input") {
        error() << "[RealDigitizer] Unknown HitEnergyContent '" << hitEnergy
                << "'. Allowed values: 'digitized', 'input'." << endmsg;
        return StatusCode::FAILURE;
      }

      if (!setupMIPCalibration()) return StatusCode::FAILURE;
      if (!setupTriggerTables()) return StatusCode::FAILURE;

      if (mode == "real") {
        m_digitizedEnergyHandle =
            std::make_unique<k4FWCore::DataHandle<podio::UserDataCollection<float>>>(
                m_digitizedEnergyName.value(), Gaudi::DataHandle::Writer, this);
        m_digitizedTimeHandle =
            std::make_unique<k4FWCore::DataHandle<podio::UserDataCollection<float>>>(
                m_digitizedTimeName.value(), Gaudi::DataHandle::Writer, this);
        m_digitizedFastHandle =
            std::make_unique<k4FWCore::DataHandle<podio::UserDataCollection<float>>>(
                m_digitizedFastName.value(), Gaudi::DataHandle::Writer, this);
        m_digitizedTriggerHandle =
            std::make_unique<k4FWCore::DataHandle<podio::UserDataCollection<int>>>(
                m_digitizedTriggerName.value(), Gaudi::DataHandle::Writer, this);
        if (m_hitSelection.value() != "cell" && m_hitSelection.value() != "chip") {
          error() << "[RealDigitizer] HitSelection must be 'cell' or 'chip', got '"
                  << m_hitSelection.value() << "'" << endmsg;
          return StatusCode::FAILURE;
        }
        info() << "[RealDigitizer] DigitizationMode='real': using contribution-level "
                  "cell shaping fast_search." << endmsg;
        info() << "[RealDigitizer] HitSelection=" << m_hitSelection.value()
               << (m_hitSelection.value() == "chip"
                       ? " (every sampled cell is kept with its trigger flag; the chip-level "
                         "decision is AdcDigitizer's)"
                       : " (only cells whose own discriminator fired)")
               << endmsg;
        info() << "[RealDigitizer] HitEnergyContent=" << hitEnergy
               << "  DigitizedEnergyCollection=" << m_digitizedEnergyName.value()
               << "  DigitizedTimeCollection=" << m_digitizedTimeName.value()
               << endmsg;
      }
      info() << "[RealDigitizer] Mode=" << mode
             << "  InputEnergyUnit=" << unit
             << "  Threshold=" << m_threshold.value() << " MIP"
             << "  FastNoise=" << m_fastNoiseMIP.value() << " MIP"
             << "  TriggerEfficiency=" << m_triggerEfficiency.value() << endmsg;
      return Gaudi::Algorithm::initialize();
    } catch (const std::exception& e) {
      error() << "[RealDigitizer] Exception in initialize(): " << e.what() << endmsg;
      return StatusCode::FAILURE;
    } catch (...) {
      error() << "[RealDigitizer] Unknown exception in initialize()." << endmsg;
      return StatusCode::FAILURE;
    }
  }

  StatusCode execute(const EventContext&) const override {
    try {
      const auto* input  = m_inputHandle->get();
      auto*       output = m_outputHandle->createAndPut();

      const long long evtNum  = m_eventCount.fetch_add(1);
      const bool      doPrint = (evtNum % m_debugFreq == 0);
      const std::string& mode = m_mode.value();
      const bool inputIsGeV = (m_inputEnergyUnit.value() == "GeV");

      int n_pass = 0;
      if (mode == "simple") {
        // Simple MIP threshold — same as BasicDigitizer.
        for (const auto& hit : *input) {
          const double mipValue = mipValueFor(hit.getCellID());
          const double invMip   = inputIsGeV ? (1.0 / mipValue) : 1.0;
          const float  energyMip = static_cast<float>(hit.getEnergy() * invMip);
          if (energyMip > m_threshold) {
            auto nh = output->create();
            nh.setCellID(hit.getCellID());
            nh.setEnergy(energyMip);
            nh.setPosition(hit.getPosition());
            for (const auto& contrib : hit.getContributions()) {
              nh.addToContributions(contrib);
            }
            ++n_pass;

            if (doPrint) {
              const auto& pos = hit.getPosition();
              debug() << "  Hit: energy=" << hit.getEnergy() << " " << m_inputEnergyUnit.value()
                      << "  => " << energyMip << " MIP"
                      << "  contributions=" << hit.contributions_size()
                      << "  pos=(" << pos.x << ", " << pos.y << ", " << pos.z << ") mm"
                      << endmsg;
            }
          }
        }
      } else if (mode == "real") {
        auto* digitizedEnergy = m_digitizedEnergyHandle->createAndPut();
        auto* digitizedTime = m_digitizedTimeHandle->createAndPut();
        auto* digitizedFast = m_digitizedFastHandle->createAndPut();
        auto* digitizedTrigger = m_digitizedTriggerHandle->createAndPut();
        const bool writeDigitizedToHit = (m_hitEnergyContent.value() == "digitized");
        const bool chipSelection = (m_hitSelection.value() == "chip");

        siwecal::CellShapingConfig cfg;
        cfg.mipThreshold = m_threshold.value();
        cfg.delayNs = m_delayNs.value();
        cfg.tauFastNs = m_tauFastNs.value();
        cfg.tauSlowNs = m_tauSlowNs.value();
        cfg.orderFast = m_orderFast.value();
        cfg.orderSlow = m_orderSlow.value();
        cfg.fastWindowNs = m_fastWindowNs.value();
        cfg.slowWindowNs = m_slowWindowNs.value();
        cfg.fastNoiseMIP = m_fastNoiseMIP.value();
        cfg.slowNoiseMIP = m_slowNoiseMIP.value();
        cfg.peakSearchBins = m_peakSearchBins.value();
        cfg.refineIterations = m_refineIterations.value();
        cfg.triggerSearchBins = m_triggerSearchBins.value();
        cfg.triggerEfficiency = m_triggerEfficiency.value();

        std::size_t hitIndex = 0;
        for (const auto& hit : *input) {
          // The shaping works in GeV steps normalised by this cell's MIP value,
          // and the trigger model is this cell's layer's (per-layer vectors, or
          // the set's scalars when they are empty).
          cfg.mipValueGeV = mipValueFor(hit.getCellID());
          setTriggerModel(hit, cfg);
          const double hitToGeV = inputIsGeV ? 1.0 : cfg.mipValueGeV;

          std::vector<double> stepEnergyGeV;
          std::vector<double> stepTimeNs;
          stepEnergyGeV.reserve(hit.contributions_size());
          stepTimeNs.reserve(hit.contributions_size());

          for (const auto& contrib : hit.getContributions()) {
            stepEnergyGeV.push_back(contrib.getEnergy());
            stepTimeNs.push_back(contrib.getTime() - m_timeOffsetNs.value());
          }

          if (stepEnergyGeV.empty()) {
            stepEnergyGeV.push_back(hit.getEnergy() * hitToGeV);
            stepTimeNs.push_back(0.0);
          }

          const auto seed = static_cast<std::uint64_t>(m_randomSeed.value())
              + 0x9e3779b97f4a7c15ULL * static_cast<std::uint64_t>(evtNum + 1)
              + static_cast<std::uint64_t>(hitIndex);
          std::mt19937_64 rng(seed);
          const auto fastSearch =
              siwecal::fastSearchCellSteps(stepEnergyGeV, stepTimeNs, cfg, rng);
          // The chip's hit selection is the fast discriminator alone: once it
          // fires, the slow sample is recorded at hold whatever its value, and
          // the data keep it too (`bestScaPerChannel` asks for hitbit_high, not
          // for an ADC level).  Cutting the slow sample at the threshold as well
          // multiplied the measured S-curve by a second, sharper one, and moved
          // the simulated 50% point from threshold_adc/adc_per_mip (0.94 MIP)
          // up to 1.0 MIP with 2/3 of the width -- exactly at the muon peak,
          // where it costs the most hits.  Only a non-positive amplitude is
          // dropped: it would quantise to nothing and carries no information.
          const bool triggered = fastSearch.triggerTime >= 0.0;
          // HitSelection='chip': a cell that did not fire is still sampled when
          // the chip's hold lands (the hold is common to the chip), so it is
          // written with its flag and AdcDigitizer decides, chip by chip, whether
          // it becomes a hit -- the test beam's EcalEventBuilder does the same
          // with HitSelection='adc'.  Cells below the shaping's pre-cut (far
          // under the threshold) carry no sample and are dropped either way.
          const bool passThreshold =
              (triggered || chipSelection) && fastSearch.slowSignalSample > 0.0;

          if (passThreshold) {
            auto nh = output->create();
            nh.setCellID(hit.getCellID());
            nh.setEnergy(writeDigitizedToHit
                             ? static_cast<float>(fastSearch.slowSignalSample)
                             : hit.getEnergy());
            nh.setPosition(hit.getPosition());
            for (const auto& contrib : hit.getContributions()) {
              nh.addToContributions(contrib);
            }
            digitizedEnergy->create() = static_cast<float>(fastSearch.slowSignalSample);
            digitizedTime->create() = static_cast<float>(fastSearch.triggerTime);
            digitizedFast->create() = static_cast<float>(fastSearch.fastPeakSample);
            digitizedTrigger->create() = triggered ? 1 : 0;
            ++n_pass;
          }

          if (doPrint) {
            const auto& pos = hit.getPosition();
            const double rawEnergyGeV =
                std::accumulate(stepEnergyGeV.begin(), stepEnergyGeV.end(), 0.0);
            debug() << "  Hit: energy=" << hit.getEnergy() << " " << m_inputEnergyUnit.value()
                    << "  mip=" << cfg.mipValueGeV << " GeV/MIP"
                    << "  raw=" << rawEnergyGeV / cfg.mipValueGeV << " MIP"
                    << "  steps=" << stepEnergyGeV.size()
                    << "  trigger=" << fastSearch.triggerTime << " ns"
                    << "  slow_sample=" << fastSearch.slowSignalSample << " MIP"
                    << "  pass=" << passThreshold
                    << "  pos=(" << pos.x << ", " << pos.y << ", " << pos.z << ") mm"
                    << endmsg;
          }
          ++hitIndex;
        }
      }

      debug() << "RealDigitizer [" << mode << "]: "
              << input->size() << " in, " << n_pass << " passing threshold" << endmsg;
      return StatusCode::SUCCESS;
    } catch (const std::exception& e) {
      error() << "[RealDigitizer] Exception in execute(): " << e.what() << endmsg;
      return StatusCode::FAILURE;
    } catch (...) {
      error() << "[RealDigitizer] Unknown exception in execute()." << endmsg;
      return StatusCode::FAILURE;
    }
  }

  StatusCode finalize() override {
    try {
      m_inputHandle.reset();
      m_outputHandle.reset();
      m_digitizedEnergyHandle.reset();
      m_digitizedTimeHandle.reset();
      m_digitizedFastHandle.reset();
      m_digitizedTriggerHandle.reset();
      m_decoder.reset();
      if (m_useTriggerTables) {
        info() << "[RealDigitizer] Trigger model used per cell: channel table " << m_nTrigChannel.load()
               << ", chip table " << m_nTrigChip.load() << ", per-layer fallback " << m_nTrigLayer.load()
               << "; with the cell's own gain " << m_nTrigGain.load() << endmsg;
      }
      return Gaudi::Algorithm::finalize();
    } catch (const std::exception& e) {
      error() << "[RealDigitizer] Exception in finalize(): " << e.what() << endmsg;
      return StatusCode::FAILURE;
    } catch (...) {
      error() << "[RealDigitizer] Unknown exception in finalize()." << endmsg;
      return StatusCode::FAILURE;
    }
  }

private:
  // --------------------------------------------------------------------- //
  // MIP calibration — same contract as GeV2MIPConversion, so a pipeline can
  // hand both algorithms the identical MIPValues list.
  // --------------------------------------------------------------------- //
  bool setupMIPCalibration() {
    const auto& mipVec = m_MIPValues.value();
    if (!mipVec.empty()) {
      const int nLayers = m_nLayers.value();
      if (nLayers > 0 && static_cast<int>(mipVec.size()) != nLayers) {
        error() << "[RealDigitizer] MIPValues has " << mipVec.size()
                << " entries but NLayers=" << nLayers
                << ". Number of MIP values must match number of detector layers." << endmsg;
        return false;
      }
      if (m_bitField.value().empty()) {
        error() << "[RealDigitizer] Per-layer MIPValues requires BitField "
                   "to decode the 'layer' field from CellID." << endmsg;
        return false;
      }
      for (std::size_t i = 0; i < mipVec.size(); ++i) {
        if (mipVec[i] <= 0.0) {
          error() << "[RealDigitizer] MIPValues[" << i << "]=" << mipVec[i]
                  << " is not positive." << endmsg;
          return false;
        }
      }
      m_decoder = std::make_unique<dd4hep::DDSegmentation::BitFieldCoder>(m_bitField.value());
      info() << "[RealDigitizer] Per-layer MIP mode: " << mipVec.size() << " layers." << endmsg;
      // The per-layer trigger model rides on the same layer decoding.
      for (const auto* vec : {&m_thresholdPerLayer.value(), &m_fastNoisePerLayer.value(),
                              &m_triggerEfficiencyPerLayer.value()}) {
        if (!vec->empty() && vec->size() != mipVec.size()) {
          error() << "[RealDigitizer] a per-layer trigger vector has " << vec->size()
                  << " entries but MIPValues has " << mipVec.size() << endmsg;
          return false;
        }
      }
      if (!m_thresholdPerLayer.value().empty()) {
        info() << "[RealDigitizer] Per-layer trigger model: threshold [MIP] =";
        for (double v : m_thresholdPerLayer.value()) info() << " " << v;
        info() << endmsg;
      }
      return true;
    }
    if (!m_thresholdPerLayer.value().empty() || !m_fastNoisePerLayer.value().empty() ||
        !m_triggerEfficiencyPerLayer.value().empty()) {
      error() << "[RealDigitizer] per-layer trigger vectors need per-layer MIPValues "
                 "(the layer is decoded the same way)." << endmsg;
      return false;
    }
    if (m_MIPValue.value() <= 0.0) {
      error() << "[RealDigitizer] MIPValue must be positive." << endmsg;
      return false;
    }
    info() << "[RealDigitizer] Scalar MIP mode: " << m_MIPValue.value()
           << " GeV/MIP." << endmsg;
    return true;
  }


  // --------------------------------------------------------------------- //
  // Trigger model per cell.  The turn-on of every channel and chip comes from
  // the muon calibration (siwecal-tb2026 PedestalMipCalibrator Mode=FitTrigger:
  // trigger/<th>/Trigger_{channel,chip}_highgain.txt, ADC above the pedestal).
  // The cell's chip and channel are found exactly as ChannelMapper finds them
  // (nearest pad in the slab's pad map, PadMap.hh).  Order: the channel's row,
  // else its chip's, else the per-layer vectors, else the set's scalars.
  // --------------------------------------------------------------------- //
  // In ADC above the pedestal, as measured; turned into MIP per cell with the
  // cell's own gain (AdcPerMip x gain shape), the same gain AdcDigitizer gives
  // the cell's charge, so a channel with a high gain reaches its ADC threshold
  // with less charge exactly as it does in the detector.
  struct TriggerModel {
    double muAdc, sigmaAdc, efficiency;
  };

  bool readTriggerTable(const std::string& path, bool perChannel) {
    std::ifstream fin(path);
    if (!fin) {
      error() << "[RealDigitizer] Cannot open trigger table " << path << endmsg;
      return false;
    }
    std::string line;
    std::size_t n = 0;
    while (std::getline(fin, line)) {
      if (line.empty() || line[0] == '#') continue;
      std::istringstream ss(line);
      int layer = -1, chip = -1, channel = -1;
      double mu, sigma, plateau;
      ss >> layer >> chip;
      if (perChannel) ss >> channel;
      ss >> mu >> sigma >> plateau;
      if (ss.fail() || !std::isfinite(mu) || !std::isfinite(sigma) || !std::isfinite(plateau)) continue;
      const TriggerModel t{mu, sigma, plateau};
      if (perChannel) {
        m_trigChannel[{layer, chip, channel}] = t;
      } else {
        m_trigChip[{layer, chip}] = t;
      }
      ++n;
    }
    info() << "[RealDigitizer] Trigger table " << path << ": " << n << " " << (perChannel ? "channels" : "chips")
           << endmsg;
    return n > 0;
  }

  bool setupTriggerTables() {
    const bool haveChip = !m_triggerChipTable.value().empty();
    const bool haveChannel = !m_triggerChannelTable.value().empty();
    if (!haveChip && !haveChannel) return true;
    if (!(m_triggerAdcPerMip.value() > 0.0)) {
      error() << "[RealDigitizer] TriggerAdcPerMip must be set (> 0) with the trigger tables" << endmsg;
      return false;
    }
    if (m_bitField.value().empty()) {
      error() << "[RealDigitizer] the trigger tables need BitField to decode the layer" << endmsg;
      return false;
    }
    if (!m_decoder) m_decoder = std::make_unique<dd4hep::DDSegmentation::BitFieldCoder>(m_bitField.value());
    std::string err;
    if (!siwecal::loadPadMap(m_padMapFile.value(), m_padDefault, err) ||
        !siwecal::loadPadMapOverrides(m_padMapSlabOverrides.value(), m_padOverrides, err)) {
      error() << "[RealDigitizer] " << err << endmsg;
      return false;
    }
    if (!m_triggerGainShapeTable.value().empty()) {
      std::string log;
      if (!siwecal::readMipTable(m_triggerGainShapeTable.value(), m_gainShape, err) ||
          !siwecal::normaliseGainShape(m_gainShape, m_gainShapeNormalisation.value(), m_gainShapeShrink.value(),
                                       m_gainShapeTrueSpread.value(), err, log)) {
        error() << "[RealDigitizer] trigger gain shape: " << err << endmsg;
        return false;
      }
      m_haveGainShape = true;
      info() << "[RealDigitizer] trigger thresholds in MIP per cell with the gain shape of "
             << m_triggerGainShapeTable.value() << " (" << m_gainShape.nCalibrated << " channels)" << endmsg;
    }
    if (haveChip && !readTriggerTable(m_triggerChipTable.value(), false)) return false;
    if (haveChannel && !readTriggerTable(m_triggerChannelTable.value(), true)) return false;
    m_useTriggerTables = true;
    return true;
  }

  void setTriggerModel(const edm4hep::SimCalorimeterHit& hit, siwecal::CellShapingConfig& cfg) const {
    const std::uint64_t cellID = hit.getCellID();
    if (m_useTriggerTables) {
      const int slab = static_cast<int>(m_decoder->get(cellID, "layer"));
      const auto& pos = hit.getPosition();
      int chip = 0, channel = 0;
      float dist2 = 0.f;
      siwecal::nearestPad(m_padDefault, m_padOverrides, slab, pos.x, pos.y, chip, channel, dist2);
      if (std::sqrt(dist2) <= m_posTolerance.value()) {
        const TriggerModel* t = nullptr;
        if (auto it = m_trigChannel.find({slab, chip, channel}); it != m_trigChannel.end()) {
          t = &it->second;
          ++m_nTrigChannel;
        } else if (auto ic = m_trigChip.find({slab, chip}); ic != m_trigChip.end()) {
          t = &ic->second;
          ++m_nTrigChip;
        }
        if (t != nullptr) {
          double adcPerMip = m_triggerAdcPerMip.value();
          if (m_haveGainShape) {
            const double g = m_gainShape.at(slab, chip, channel);
            if (std::isfinite(g) && g > 0.0) {
              adcPerMip *= g;
              ++m_nTrigGain;
            }
          }
          cfg.mipThreshold = t->muAdc / adcPerMip;
          cfg.fastNoiseMIP = t->sigmaAdc / adcPerMip;
          cfg.triggerEfficiency = t->efficiency;
          return;
        }
      }
    }
    if (m_useTriggerTables) ++m_nTrigLayer;
    cfg.mipThreshold = perLayerOr(m_thresholdPerLayer.value(), cellID, m_threshold.value());
    cfg.fastNoiseMIP = perLayerOr(m_fastNoisePerLayer.value(), cellID, m_fastNoiseMIP.value());
    cfg.triggerEfficiency = perLayerOr(m_triggerEfficiencyPerLayer.value(), cellID, m_triggerEfficiency.value());
  }

  Gaudi::Property<std::string> m_triggerChipTable{
      this, "TriggerChipTable", "",
      "Per-chip trigger turn-on, 'layer chip mu_adc sigma_adc plateau ...' (siwecal-tb2026 "
      "PedestalMipCalibrator Mode=FitTrigger); empty = per-layer vectors / scalars"};
  Gaudi::Property<std::string> m_triggerChannelTable{
      this, "TriggerChannelTable", "",
      "Per-channel trigger turn-on, 'layer chip channel mu_adc sigma_adc plateau ...'; wins over the chip table"};
  Gaudi::Property<double> m_triggerAdcPerMip{
      this, "TriggerAdcPerMip", 0.0,
      "ADC per MIP of the threshold set, to turn the tables' ADC into the MIP the shaping works in"};
  Gaudi::Property<std::string> m_triggerGainShapeTable{
      this, "TriggerGainShapeTable", "",
      "MIP table whose gain shape turns each cell's ADC threshold into MIP (AdcDigitizer's GainShapeThreshold "
      "table); empty = TriggerAdcPerMip for every cell"};
  Gaudi::Property<std::string> m_gainShapeNormalisation{
      this, "GainShapeNormalisation", "per-slab", "As AdcDigitizer.GainShapeNormalisation"};
  Gaudi::Property<bool> m_gainShapeShrink{this, "GainShapeShrink", false, "As AdcDigitizer.GainShapeShrink"};
  Gaudi::Property<double> m_gainShapeTrueSpread{
      this, "GainShapeTrueSpread", 0.0, "As AdcDigitizer.GainShapeTrueSpread"};
  Gaudi::Property<std::string> m_padMapFile{
      this, "PadMapFile", "mappings/fev10_rotate_chip_channel_x_y_mapping.txt",
      "Default pad map (chip/channel -> x,y), as ChannelMapper.PadMapFile"};
  Gaudi::Property<std::vector<std::string>> m_padMapSlabOverrides{
      this, "PadMapSlabOverrides", {}, "Per-slab pad maps 'slab:path', as ChannelMapper.PadMapSlabOverrides"};
  Gaudi::Property<double> m_posTolerance{
      this, "PositionTolerance", 4.0,
      "Max distance [mm] to the nearest pad, as ChannelMapper.PositionTolerance; beyond it the cell keeps the "
      "per-layer model"};

  bool m_useTriggerTables = false;
  bool m_haveGainShape = false;
  siwecal::MipTable m_gainShape;
  mutable std::atomic<long long> m_nTrigChannel{0}, m_nTrigChip{0}, m_nTrigLayer{0}, m_nTrigGain{0};
  std::vector<siwecal::PadEntry> m_padDefault;
  std::map<int, std::vector<siwecal::PadEntry>> m_padOverrides;
  std::map<std::tuple<int, int, int>, TriggerModel> m_trigChannel;
  std::map<std::pair<int, int>, TriggerModel> m_trigChip;

  double mipValueFor(std::uint64_t cellID) const {
    return perLayerOr(m_MIPValues.value(), cellID, m_MIPValue.value());
  }

  // A per-layer vector's entry for this cell, or `fallback` when the vector is
  // empty (scalar mode).  The trigger model uses it exactly like the MIP value:
  // slab 12 of the 2026 stack ran at its own, higher discriminator (DAC 243
  // against 215-230), and a slab that triggers later needs its own threshold and
  // width, not the set's.
  double perLayerOr(const std::vector<double>& vec, std::uint64_t cellID, double fallback) const {
    if (vec.empty()) return fallback;
    const int layer = static_cast<int>(m_decoder->get(cellID, "layer"));
    if (layer < 0 || layer >= static_cast<int>(vec.size())) {
      warning() << "[RealDigitizer] Layer " << layer
                << " out of range [0," << vec.size()
                << "), using first entry." << endmsg;
      return vec[0];
    }
    return vec[layer];
  }

  Gaudi::Property<std::string> m_inputName{
      this, "InputCollection", "SiPadHitsMIP",
      "Input SimCalorimeterHit collection (in MIP units)"};
  Gaudi::Property<std::string> m_outputName{
      this, "OutputCollection", "SiPadHitsDigi",
      "Output digitized SimCalorimeterHit collection"};
  Gaudi::Property<std::string> m_digitizedEnergyName{
      this, "DigitizedEnergyCollection", "SiPadHitsDigiDigitizedEnergy",
      "Output UserDataCollection<float> with shaped slow-sample amplitudes [MIP]"};
  Gaudi::Property<std::string> m_digitizedTimeName{
      this, "DigitizedTimeCollection", "SiPadHitsDigiDigitizedTime",
      "Output UserDataCollection<float> with fast trigger times [ns]"};
  Gaudi::Property<std::string> m_digitizedTriggerName{
      this, "DigitizedTriggerCollection", "SiPadHitsDigiDigitizedTrigger",
      "Output UserDataCollection<int>, 1 if the cell's own discriminator fired "
      "(the test beam's hitbit_high), 0 if it was only sampled"};
  Gaudi::Property<std::string> m_hitSelection{
      this, "HitSelection", "cell",
      "'cell': a cell becomes a hit only if its own discriminator fired (the "
      "historical behaviour, EcalEventBuilder HitSelection='hitbit').  'chip': "
      "every sampled cell is written with its trigger flag and AdcDigitizer keeps, "
      "on chips where at least one cell fired, the cells that fired or exceed its "
      "AdcHitThreshold (EcalEventBuilder HitSelection='adc')"};
  Gaudi::Property<std::string> m_digitizedFastName{
      this, "DigitizedFastCollection", "SiPadHitsDigiDigitizedFast",
      "Output UserDataCollection<float> with fast-channel peak amplitudes [MIP] "
      "-- what the trigger discriminator sees, to be compared with the test "
      "beam's hitbit_high turn-on"};
  Gaudi::Property<double> m_threshold{
      this, "Threshold", 0.5,
      "Minimum energy [MIP] to keep a hit"};
  Gaudi::Property<std::string> m_inputEnergyUnit{
      this, "InputEnergyUnit", "MIP",
      "Input hit energy unit: 'MIP' or 'GeV'"};
  Gaudi::Property<std::string> m_hitEnergyContent{
      this, "HitEnergyContent", "digitized",
      "What mode='real' writes into the output hit energy: 'digitized' (shaped "
      "slow-sample amplitude [MIP]) or 'input' (original input hit energy)"};
  Gaudi::Property<double> m_MIPValue{
      this, "MIPValue", 0.0002,
      "MIP calibration value [GeV/MIP] (scalar mode)"};
  Gaudi::Property<std::vector<double>> m_MIPValues{
      this, "MIPValues", {},
      "Per-layer MIP calibration [GeV/MIP]. Non-empty overrides MIPValue and "
      "requires BitField"};
  Gaudi::Property<int> m_nLayers{
      this, "NLayers", 0,
      "Expected length of MIPValues (0 skips the check)"};
  Gaudi::Property<std::vector<double>> m_thresholdPerLayer{
      this, "ThresholdPerLayer", {},
      "Per-layer discriminator level [MIP]; empty = Threshold for every layer.  Same "
      "length as MIPValues"};
  Gaudi::Property<std::vector<double>> m_fastNoisePerLayer{
      this, "FastNoiseMIPPerLayer", {},
      "Per-layer turn-on width [MIP]; empty = FastNoiseMIP for every layer"};
  Gaudi::Property<std::vector<double>> m_triggerEfficiencyPerLayer{
      this, "TriggerEfficiencyPerLayer", {},
      "Per-layer plateau efficiency; empty = TriggerEfficiency for every layer"};
  Gaudi::Property<std::string> m_bitField{
      this, "BitField", "system:8,layer:8,slice:5,x:9,y:9",
      "CellID bitfield used to decode 'layer' in per-layer MIP mode"};
  Gaudi::Property<std::string> m_mode{
      this, "DigitizationMode", "simple",
      "Digitization mode: 'simple' (MIP threshold cut) or 'real' "
      "(cell shaping fast_search)"};
  Gaudi::Property<int> m_debugFreq{
      this, "DebugFrequency", 500,
      "Print per-hit debug info every N events"};
  Gaudi::Property<double> m_delayNs{
      this, "DelayNs", 160.0,
      "Delay [ns] between fast trigger time and slow sample"};
  Gaudi::Property<double> m_tauFastNs{
      this, "TauFastNs", 30.0,
      "Fast CR-RC shaping time [ns]"};
  Gaudi::Property<double> m_tauSlowNs{
      this, "TauSlowNs", 180.0,
      "Slow CR-RC shaping time [ns]"};
  Gaudi::Property<int> m_orderFast{
      this, "OrderFast", 2,
      "Fast CR-RC shaping order"};
  Gaudi::Property<int> m_orderSlow{
      this, "OrderSlow", 2,
      "Slow CR-RC shaping order"};
  Gaudi::Property<double> m_timeOffsetNs{
      this, "TimeOffsetNs", 0.0,
      "Subtracted from every contribution time [ns]: the extra flight time when the beam is born further "
      "upstream than the -2650 mm the shaping windows were set for (job3_digitize.py: from BEAM_Z_MM)"};
  Gaudi::Property<double> m_fastWindowNs{
      this, "FastWindowNs", 200.0,
      "Fast peak search window [ns]"};
  Gaudi::Property<double> m_slowWindowNs{
      this, "SlowWindowNs", 500.0,
      "Slow peak search window [ns]"};
  Gaudi::Property<double> m_fastNoiseMIP{
      this, "FastNoiseMIP", 1.0 / 30.0,
      "Fast peak Gaussian noise sigma [MIP]"};
  Gaudi::Property<double> m_slowNoiseMIP{
      this, "SlowNoiseMIP", 1.0 / 12.0,
      "Slow peak/sample Gaussian noise sigma [MIP]"};
  Gaudi::Property<int> m_peakSearchBins{
      this, "PeakSearchBins", 64,
      "Coarse bins for fast/slow peak search"};
  Gaudi::Property<int> m_refineIterations{
      this, "RefineIterations", 48,
      "Iterations for peak/threshold refinement"};
  Gaudi::Property<int> m_triggerSearchBins{
      this, "TriggerSearchBins", 64,
      "Coarse bins for fast-threshold trigger search"};
  Gaudi::Property<double> m_triggerEfficiency{
      this, "TriggerEfficiency", 1.0,
      "Efficiency above the threshold [0,1].  The test beam's discriminator "
      "plateaus at 0.85-0.96 depending on the threshold set, not at 1; measured "
      "by analysis/trigger_turnon.py and configured per th in "
      "mappings/trigger_turnon.yml"};
  Gaudi::Property<unsigned long long> m_randomSeed{
      this, "RandomSeed", 5489ULL,
      "Base seed for shaping noise"};

  mutable std::unique_ptr<k4FWCore::DataHandle<edm4hep::SimCalorimeterHitCollection>> m_inputHandle;
  mutable std::unique_ptr<k4FWCore::DataHandle<edm4hep::SimCalorimeterHitCollection>> m_outputHandle;
  mutable std::unique_ptr<k4FWCore::DataHandle<podio::UserDataCollection<float>>> m_digitizedEnergyHandle;
  mutable std::unique_ptr<k4FWCore::DataHandle<podio::UserDataCollection<float>>> m_digitizedTimeHandle;
  mutable std::unique_ptr<k4FWCore::DataHandle<podio::UserDataCollection<float>>> m_digitizedFastHandle;
  mutable std::unique_ptr<k4FWCore::DataHandle<podio::UserDataCollection<int>>> m_digitizedTriggerHandle;
  std::unique_ptr<dd4hep::DDSegmentation::BitFieldCoder> m_decoder;
  mutable std::atomic<long long> m_eventCount{0};
};

DECLARE_COMPONENT(RealDigitizer)

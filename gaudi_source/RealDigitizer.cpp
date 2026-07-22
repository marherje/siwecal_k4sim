#include "Gaudi/Algorithm.h"
#include "GaudiKernel/MsgStream.h"
#include "k4FWCore/DataHandle.h"
#include "CellShaping.hh"
#include "edm4hep/MCParticle.h"
#include "edm4hep/SimCalorimeterHitCollection.h"
#include "podio/UserDataCollection.h"

#include <atomic>
#include <cstdint>
#include <memory>
#include <numeric>
#include <random>
#include <vector>

// Realistic SiPad digitizer — placeholder implementation.
//
// DigitizationMode:
//   "simple" — identical to BasicDigitizer: apply an energy threshold cut.
//   "real"   — run contribution-level CR-RC fast_search shaping without
//               building a full waveform array.
//
// InputEnergyUnit:
//   "MIP" — input hit energy is already in MIP units.
//   "GeV" — input hit energy is converted with MIPValue [GeV/MIP].
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
      if (m_MIPValue.value() <= 0.0) {
        error() << "[RealDigitizer] MIPValue must be positive." << endmsg;
        return StatusCode::FAILURE;
      }
      if (mode == "real") {
        m_digitizedEnergyHandle =
            std::make_unique<k4FWCore::DataHandle<podio::UserDataCollection<float>>>(
                m_digitizedEnergyName.value(), Gaudi::DataHandle::Writer, this);
        m_digitizedTimeHandle =
            std::make_unique<k4FWCore::DataHandle<podio::UserDataCollection<float>>>(
                m_digitizedTimeName.value(), Gaudi::DataHandle::Writer, this);
        info() << "[RealDigitizer] DigitizationMode='real': using contribution-level "
                  "cell shaping fast_search." << endmsg;
        info() << "[RealDigitizer] DigitizedEnergyCollection="
               << m_digitizedEnergyName.value()
               << "  DigitizedTimeCollection=" << m_digitizedTimeName.value()
               << endmsg;
      }
      info() << "[RealDigitizer] Mode=" << mode
             << "  InputEnergyUnit=" << unit
             << "  MIPValue=" << m_MIPValue.value() << " GeV/MIP"
             << "  Threshold=" << m_threshold.value() << " MIP" << endmsg;
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
      const double invMip = inputIsGeV ? (1.0 / m_MIPValue.value()) : 1.0;
      const double hitToGeV = inputIsGeV ? 1.0 : m_MIPValue.value();

      int n_pass = 0;
      if (mode == "simple") {
        // Simple MIP threshold — same as BasicDigitizer.
        for (const auto& hit : *input) {
          const float energyMip = static_cast<float>(hit.getEnergy() * invMip);
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

              int contribIdx = 0;
              for (const auto& contrib : hit.getContributions()) {
                const auto particle = contrib.getParticle();
                const int pdg = particle.isAvailable() ? particle.getPDG() : contrib.getPDG();
                const auto& step = contrib.getStepPosition();
                debug() << "    contrib=" << contribIdx
                        << " energy=" << contrib.getEnergy() << " GeV"
                        << " => " << (contrib.getEnergy() / m_MIPValue.value()) << " MIP"
                        << " time=" << contrib.getTime() << " ns"
                        << " pdg=" << pdg
                        << " step=(" << step.x << ", " << step.y << ", " << step.z << ") mm"
                        << endmsg;
                ++contribIdx;
              }
            }
          }
        }
      } else if (mode == "real") {
        auto* digitizedEnergy = m_digitizedEnergyHandle->createAndPut();
        auto* digitizedTime = m_digitizedTimeHandle->createAndPut();

        siwecal::CellShapingConfig cfg;
        cfg.mipValueGeV = m_MIPValue.value();
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

        std::size_t hitIndex = 0;
        for (const auto& hit : *input) {
          std::vector<double> stepEnergyGeV;
          std::vector<double> stepTimeNs;
          stepEnergyGeV.reserve(hit.contributions_size());
          stepTimeNs.reserve(hit.contributions_size());

          for (const auto& contrib : hit.getContributions()) {
            stepEnergyGeV.push_back(contrib.getEnergy());
            stepTimeNs.push_back(contrib.getTime());
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
          const bool passThreshold =
              fastSearch.triggerTime >= 0.0 &&
              fastSearch.slowSignalSample >= cfg.mipThreshold;

          if (passThreshold) {
            auto nh = output->create();
            nh.setCellID(hit.getCellID());
            nh.setEnergy(hit.getEnergy());
            nh.setPosition(hit.getPosition());
            for (const auto& contrib : hit.getContributions()) {
              nh.addToContributions(contrib);
            }
            digitizedEnergy->create() = static_cast<float>(fastSearch.slowSignalSample);
            digitizedTime->create() = static_cast<float>(fastSearch.triggerTime);
            ++n_pass;
          }

          if (doPrint) {
            const auto& pos = hit.getPosition();
            const double rawEnergyGeV =
                std::accumulate(stepEnergyGeV.begin(), stepEnergyGeV.end(), 0.0);
            debug() << "  Hit: energy=" << hit.getEnergy() << " " << m_inputEnergyUnit.value()
                    << "  fast_search raw=" << rawEnergyGeV / m_MIPValue.value()
                    << " MIP"
                    << "  steps=" << stepEnergyGeV.size()
                    << "  fast_search trigger=" << fastSearch.triggerTime << " ns"
                    << "  digitized_time=" << fastSearch.triggerTime << " ns"
                    << "  fast_search slow_sample=" << fastSearch.slowSignalSample
                    << " MIP"
                    << "  digitized_energy=" << fastSearch.slowSignalSample << " MIP"
                    << "  fast_search pass=" << passThreshold
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
  Gaudi::Property<double> m_threshold{
      this, "Threshold", 0.5,
      "Minimum energy [MIP] to keep a hit"};
  Gaudi::Property<std::string> m_inputEnergyUnit{
      this, "InputEnergyUnit", "MIP",
      "Input hit energy unit: 'MIP' or 'GeV'"};
  Gaudi::Property<double> m_MIPValue{
      this, "MIPValue", 0.0002,
      "MIP calibration value [GeV/MIP] used when InputEnergyUnit='GeV'"};
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
  Gaudi::Property<unsigned long long> m_randomSeed{
      this, "RandomSeed", 5489ULL,
      "Base seed for shaping noise"};

  mutable std::unique_ptr<k4FWCore::DataHandle<edm4hep::SimCalorimeterHitCollection>> m_inputHandle;
  mutable std::unique_ptr<k4FWCore::DataHandle<edm4hep::SimCalorimeterHitCollection>> m_outputHandle;
  mutable std::unique_ptr<k4FWCore::DataHandle<podio::UserDataCollection<float>>> m_digitizedEnergyHandle;
  mutable std::unique_ptr<k4FWCore::DataHandle<podio::UserDataCollection<float>>> m_digitizedTimeHandle;
  mutable std::atomic<long long> m_eventCount{0};
};

DECLARE_COMPONENT(RealDigitizer)

#include "Gaudi/Algorithm.h"
#include "GaudiKernel/MsgStream.h"
#include "k4FWCore/DataHandle.h"
#include "edm4hep/SimCalorimeterHitCollection.h"

#include <atomic>
#include <memory>

// Realistic SiPad digitizer — placeholder implementation.
//
// DigitizationMode:
//   "simple" — identical to BasicDigitizer: apply an energy threshold cut.
//   "real"   — placeholder for future realistic digitization (charge sharing,
//               noise, ToT response, etc.). Currently falls back to "simple".
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
      if (mode == "real") {
        warning() << "[RealDigitizer] DigitizationMode='real' is a placeholder. "
                     "Running simple MIP threshold cut." << endmsg;
      }
      info() << "[RealDigitizer] Mode=" << mode
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

      int n_pass = 0;
      if (mode == "simple" || mode == "real") {
        // Simple MIP threshold — same as BasicDigitizer.
        // TODO (mode=="real"): add charge sharing, noise, ToT response, etc.
        for (const auto& hit : *input) {
          if (hit.getEnergy() > m_threshold) {
            auto nh = output->create();
            nh.setCellID(hit.getCellID());
            nh.setEnergy(hit.getEnergy());
            nh.setPosition(hit.getPosition());
            ++n_pass;

            if (doPrint) {
              const auto& pos = hit.getPosition();
              debug() << "  Hit: energy=" << hit.getEnergy() << " MIP"
                      << "  pos=(" << pos.x << ", " << pos.y << ", " << pos.z << ") mm"
                      << endmsg;
            }
          }
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
  Gaudi::Property<double> m_threshold{
      this, "Threshold", 0.5,
      "Minimum energy [MIP] to keep a hit"};
  Gaudi::Property<std::string> m_mode{
      this, "DigitizationMode", "simple",
      "Digitization mode: 'simple' (MIP threshold cut) or 'real' (placeholder)"};
  Gaudi::Property<int> m_debugFreq{
      this, "DebugFrequency", 500,
      "Print per-hit debug info every N events"};

  mutable std::unique_ptr<k4FWCore::DataHandle<edm4hep::SimCalorimeterHitCollection>> m_inputHandle;
  mutable std::unique_ptr<k4FWCore::DataHandle<edm4hep::SimCalorimeterHitCollection>> m_outputHandle;
  mutable std::atomic<long long> m_eventCount{0};
};

DECLARE_COMPONENT(RealDigitizer)

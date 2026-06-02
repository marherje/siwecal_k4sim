#include "Gaudi/Algorithm.h"
#include "GaudiKernel/MsgStream.h"
#include "k4FWCore/DataHandle.h"
#include "edm4hep/SimCalorimeterHitCollection.h"

#include <atomic>
#include <memory>

class GeV2MIPConversion : public Gaudi::Algorithm {
public:
  GeV2MIPConversion(const std::string& name, ISvcLocator* svcLoc)
      : Gaudi::Algorithm(name, svcLoc) {}

  StatusCode initialize() override {
    try {
      StatusCode sc = Gaudi::Algorithm::initialize();
      if (sc.isFailure()) return sc;
      // DataHandles are built here, AFTER Gaudi has applied all Python-set
      // property values, so m_inputName / m_outputName hold the user values.
      m_inputHandle  = std::make_unique<k4FWCore::DataHandle<edm4hep::SimCalorimeterHitCollection>>(
          m_inputName.value(),  Gaudi::DataHandle::Reader, this);
      m_outputHandle = std::make_unique<k4FWCore::DataHandle<edm4hep::SimCalorimeterHitCollection>>(
          m_outputName.value(), Gaudi::DataHandle::Writer, this);
      m_invMIP = 1.0 / m_MIPValue;
      return sc;
    } catch (const std::exception& e) {
      error() << "[GeV2MIPConversion] Exception in initialize(): " << e.what() << endmsg;
      return StatusCode::FAILURE;
    } catch (...) {
      error() << "[GeV2MIPConversion] Unknown exception in initialize()." << endmsg;
      return StatusCode::FAILURE;
    }
  }

  StatusCode execute(const EventContext&) const override {
    try {
      const auto* input  = m_inputHandle->get();
      auto*       output = m_outputHandle->createAndPut();

      const long long evtNum = m_eventCount.fetch_add(1);
      const bool doPrint = (evtNum % m_debugFreq == 0);

      for (const auto& hit : *input) {
        auto nh = output->create();
        nh.setCellID(hit.getCellID());
        const auto& pos = hit.getPosition();
        const float mipEnergy = static_cast<float>(hit.getEnergy() * m_invMIP);
        nh.setEnergy(mipEnergy);
        nh.setPosition(pos);

        if (doPrint) {
          debug() << "  Hit: energy=" << hit.getEnergy() << " GeV"
                  << "  pos=(" << pos.x << ", " << pos.y << ", " << pos.z << ") mm"
                  << "  => " << mipEnergy << " MIP"
                  << endmsg;
        }
      }

      if (doPrint) {
        info() << "GeV2MIPConversion [evt " << evtNum << "]: "
               << input->size() << " hits processed (MIPValue=" << m_MIPValue << " GeV)"
               << endmsg;
      }

      return StatusCode::SUCCESS;
    } catch (const std::exception& e) {
      error() << "[GeV2MIPConversion] Exception in execute(): " << e.what() << endmsg;
      return StatusCode::FAILURE;
    } catch (...) {
      error() << "[GeV2MIPConversion] Unknown exception in execute()." << endmsg;
      return StatusCode::FAILURE;
    }
  }

  StatusCode finalize() override {
    try {
      m_inputHandle.reset();
      m_outputHandle.reset();
      return Gaudi::Algorithm::finalize();
    } catch (const std::exception& e) {
      error() << "[GeV2MIPConversion] Exception in finalize(): " << e.what() << endmsg;
      return StatusCode::FAILURE;
    } catch (...) {
      error() << "[GeV2MIPConversion] Unknown exception in finalize()." << endmsg;
      return StatusCode::FAILURE;
    }
  }

private:
  // Python-configurable properties
  Gaudi::Property<std::string> m_inputName{
      this, "InputCollection", "SiPadHitsWindowed", "Input SimCalorimeterHit collection"};
  Gaudi::Property<std::string> m_outputName{
      this, "OutputCollection", "SiPadHitsMIP", "Output SimCalorimeterHit collection"};
  Gaudi::Property<double> m_MIPValue{
      this, "MIPValue", 0.000009, "MIP value for conversion (in GeV)"}; // Default is 9 keV, for 0.3mm silicon

  Gaudi::Property<int> m_debugFreq{
      this, "DebugFrequency", 500, "Print per-hit debug info every N events"};

  double m_invMIP{0.0}; // cached reciprocal of m_MIPValue, set in initialize()

  // DataHandles constructed in initialize() from the properties above
  mutable std::unique_ptr<k4FWCore::DataHandle<edm4hep::SimCalorimeterHitCollection>> m_inputHandle;
  mutable std::unique_ptr<k4FWCore::DataHandle<edm4hep::SimCalorimeterHitCollection>> m_outputHandle;
  mutable std::atomic<long long> m_eventCount{0};
};

DECLARE_COMPONENT(GeV2MIPConversion)

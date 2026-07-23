#include "Gaudi/Algorithm.h"
#include "GaudiKernel/MsgStream.h"
#include "k4FWCore/DataHandle.h"
#include "edm4hep/SimCalorimeterHitCollection.h"
#include "DDSegmentation/BitFieldCoder.h"

#include <atomic>
#include <memory>
#include <vector>

// Converts hit energies from GeV to MIP units.
// Accepts either a single scalar MIPValue (applied to all layers) or a
// per-layer vector MIPValues (one entry per detector layer).
// When MIPValues is non-empty it overrides MIPValue; NLayers is used to
// validate the vector length (set NLayers=0 to skip the check).
class GeV2MIPConversion : public Gaudi::Algorithm {
public:
  GeV2MIPConversion(const std::string& name, ISvcLocator* svcLoc)
      : Gaudi::Algorithm(name, svcLoc) {}

  StatusCode initialize() override {
    try {
      StatusCode sc = Gaudi::Algorithm::initialize();
      if (sc.isFailure()) return sc;
      m_inputHandle  = std::make_unique<k4FWCore::DataHandle<edm4hep::SimCalorimeterHitCollection>>(
          m_inputName.value(),  Gaudi::DataHandle::Reader, this);
      m_outputHandle = std::make_unique<k4FWCore::DataHandle<edm4hep::SimCalorimeterHitCollection>>(
          m_outputName.value(), Gaudi::DataHandle::Writer, this);

      const auto& mipVec = m_MIPValues.value();
      if (!mipVec.empty()) {
        // --- Per-layer mode ---
        const int nLayers = m_nLayers.value();
        if (nLayers > 0 && static_cast<int>(mipVec.size()) != nLayers) {
          error() << "[GeV2MIPConversion] MIPValues has " << mipVec.size()
                  << " entries but NLayers=" << nLayers
                  << ". Number of MIP values must match number of detector layers." << endmsg;
          return StatusCode::FAILURE;
        }
        if (m_bitField.value().empty()) {
          error() << "[GeV2MIPConversion] Per-layer MIPValues requires BitField "
                     "to decode the 'layer' field from CellID." << endmsg;
          return StatusCode::FAILURE;
        }
        m_decoder = std::make_unique<dd4hep::DDSegmentation::BitFieldCoder>(m_bitField.value());
        m_invMIPVec.resize(mipVec.size());
        for (std::size_t i = 0; i < mipVec.size(); i++) {
          if (mipVec[i] <= 0.0) {
            error() << "[GeV2MIPConversion] MIPValues[" << i << "]=" << mipVec[i]
                    << " is not positive." << endmsg;
            return StatusCode::FAILURE;
          }
          m_invMIPVec[i] = 1.0 / mipVec[i];
        }
        info() << "[GeV2MIPConversion] Per-layer MIP mode: " << mipVec.size() << " layers." << endmsg;
      } else {
        // --- Scalar mode ---
        if (m_MIPValue.value() <= 0.0) {
          error() << "[GeV2MIPConversion] MIPValue must be positive." << endmsg;
          return StatusCode::FAILURE;
        }
        m_invMIP = 1.0 / m_MIPValue.value();
        info() << "[GeV2MIPConversion] Scalar MIP mode: " << m_MIPValue.value() << " GeV/MIP." << endmsg;
      }
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

      const long long evtNum  = m_eventCount.fetch_add(1);
      const bool      doPrint = (evtNum % m_debugFreq == 0);
      const bool      perLayer = !m_invMIPVec.empty();

      for (const auto& hit : *input) {
        auto nh = output->create();
        nh.setCellID(hit.getCellID());
        nh.setPosition(hit.getPosition());

        double invMip;
        if (perLayer) {
          const int layer = static_cast<int>(m_decoder->get(hit.getCellID(), "layer"));
          if (layer < 0 || layer >= static_cast<int>(m_invMIPVec.size())) {
            warning() << "[GeV2MIPConversion] Layer " << layer
                      << " out of range [0," << m_invMIPVec.size()
                      << "), using first entry." << endmsg;
            invMip = m_invMIPVec[0];
          } else {
            invMip = m_invMIPVec[layer];
          }
        } else {
          invMip = m_invMIP;
        }

        const float mipEnergy = static_cast<float>(hit.getEnergy() * invMip);
        nh.setEnergy(mipEnergy);

        if (doPrint) {
          const auto& pos = hit.getPosition();
          debug() << "  Hit: energy=" << hit.getEnergy() << " GeV"
                  << "  pos=(" << pos.x << ", " << pos.y << ", " << pos.z << ") mm"
                  << "  => " << mipEnergy << " MIP"
                  << endmsg;
        }
      }

      if (doPrint) {
        info() << "GeV2MIPConversion [evt " << evtNum << "]: "
               << input->size() << " hits processed" << endmsg;
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
      m_decoder.reset();
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
  Gaudi::Property<std::string> m_inputName{
      this, "InputCollection", "SiPadHitsWindowed",
      "Input SimCalorimeterHit collection (in GeV)"};
  Gaudi::Property<std::string> m_outputName{
      this, "OutputCollection", "SiPadHitsMIP",
      "Output SimCalorimeterHit collection (in MIP)"};
  Gaudi::Property<double> m_MIPValue{
      this, "MIPValue", 0.000009,
      "MIP value [GeV] applied to all layers (used when MIPValues is empty)"};
  Gaudi::Property<std::vector<double>> m_MIPValues{
      this, "MIPValues", {},
      "Per-layer MIP values [GeV]. If non-empty, overrides MIPValue. "
      "Size must equal NLayers (when NLayers > 0)."};
  Gaudi::Property<int> m_nLayers{
      this, "NLayers", 0,
      "Expected number of detector layers for validation (0 = skip check)"};
  Gaudi::Property<std::string> m_bitField{
      this, "BitField", "system:8,layer:8,slice:5,x:9,y:9",
      "DD4hep BitField encoding string used to decode the 'layer' field from CellID"};
  Gaudi::Property<int> m_debugFreq{
      this, "DebugFrequency", 500, "Print per-hit debug info every N events"};

  // Scalar mode: cached 1/MIPValue
  double m_invMIP{0.0};
  // Per-layer mode: cached 1/MIPValues[layer]
  std::vector<double> m_invMIPVec;

  mutable std::unique_ptr<dd4hep::DDSegmentation::BitFieldCoder> m_decoder;
  mutable std::unique_ptr<k4FWCore::DataHandle<edm4hep::SimCalorimeterHitCollection>> m_inputHandle;
  mutable std::unique_ptr<k4FWCore::DataHandle<edm4hep::SimCalorimeterHitCollection>> m_outputHandle;
  mutable std::atomic<long long> m_eventCount{0};
};

DECLARE_COMPONENT(GeV2MIPConversion)

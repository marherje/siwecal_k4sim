/*
 * DetectorFlipper: rewrite the z coordinate of each CalorimeterHit according
 * to a configurable per-layer z-position table.
 *
 * Motivation
 * ----------
 * In the real TB2026 setup the detector is physically flipped with respect to
 * the simulation convention: slab 0 faces the back instead of the front.
 * Raw reconstructed hits therefore carry z values that are reversed relative
 * to what the shower-variable code (EcalPidTransformer) and the event viewer
 * expect.  This algorithm provides a single place to correct that: read the
 * layer index from the CellID, look it up in ZPositions, and overwrite z.
 * The result is a new hit collection in which z always increases from the
 * first active layer to the last, regardless of the physical orientation of
 * the detector.
 *
 * Usage in the pipeline
 * ---------------------
 * Add after BasicDigitizer / RealDigitizer in the job3 steering file:
 *
 *   from Configurables import DetectorFlipper
 *   flip = DetectorFlipper("DetectorFlipper_SiPad")
 *   flip.InputCollection  = "SiPadHitsDigi"
 *   flip.OutputCollection = "SiPadHitsFlipped"
 *   flip.ZPositions = [z0, z1, ..., z14]   # mm, one value per layer
 *
 * For simulation the defaults are the DD4hep z values (no physical change;
 * useful to verify the algorithm and keep a single-format output).
 * For real TB data, set ZPositions to the desired canonical z convention.
 *
 * Properties
 * ----------
 *   InputCollection   CalorimeterHit collection to read (default: SiPadHitsDigi)
 *   OutputCollection  CalorimeterHit collection to write (default: SiPadHitsFlipped)
 *   ZPositions        Per-layer z values [mm].  Must have exactly NLayers entries.
 *   BitField          DD4hep CellID encoding string (used to decode 'layer').
 *   NLayers           Expected layer count; used to validate ZPositions length
 *                     (set to 0 to skip the check).
 *   DebugFrequency    Print per-event debug info every N events.
 */
#include "Gaudi/Algorithm.h"
#include "GaudiKernel/MsgStream.h"
#include "k4FWCore/DataHandle.h"
#include "edm4hep/SimCalorimeterHitCollection.h"
#include "DDSegmentation/BitFieldCoder.h"

#include <atomic>
#include <memory>
#include <vector>

class DetectorFlipper : public Gaudi::Algorithm {
public:
  DetectorFlipper(const std::string& name, ISvcLocator* svcLoc)
      : Gaudi::Algorithm(name, svcLoc) {}

  StatusCode initialize() override {
    try {
      StatusCode sc = Gaudi::Algorithm::initialize();
      if (sc.isFailure()) return sc;

      m_inputHandle  = std::make_unique<k4FWCore::DataHandle<edm4hep::SimCalorimeterHitCollection>>(
          m_inputName.value(),  Gaudi::DataHandle::Reader, this);
      m_outputHandle = std::make_unique<k4FWCore::DataHandle<edm4hep::SimCalorimeterHitCollection>>(
          m_outputName.value(), Gaudi::DataHandle::Writer, this);

      const auto& zpos = m_zPositions.value();
      if (zpos.empty()) {
        error() << "[DetectorFlipper] ZPositions must be non-empty." << endmsg;
        return StatusCode::FAILURE;
      }
      const int nLayers = m_nLayers.value();
      if (nLayers > 0 && static_cast<int>(zpos.size()) != nLayers) {
        error() << "[DetectorFlipper] ZPositions has " << zpos.size()
                << " entries but NLayers=" << nLayers << "." << endmsg;
        return StatusCode::FAILURE;
      }

      if (m_bitField.value().empty()) {
        error() << "[DetectorFlipper] BitField must be set to decode the "
                   "'layer' field from CellID." << endmsg;
        return StatusCode::FAILURE;
      }
      m_decoder = std::make_unique<dd4hep::DDSegmentation::BitFieldCoder>(
          m_bitField.value());

      info() << "[DetectorFlipper] Configured for " << zpos.size()
             << " layers.  ZPositions[0]=" << zpos.front()
             << " mm, ZPositions[" << (zpos.size()-1) << "]=" << zpos.back()
             << " mm." << endmsg;
      return sc;
    } catch (const std::exception& e) {
      error() << "[DetectorFlipper] Exception in initialize(): " << e.what() << endmsg;
      return StatusCode::FAILURE;
    } catch (...) {
      error() << "[DetectorFlipper] Unknown exception in initialize()." << endmsg;
      return StatusCode::FAILURE;
    }
  }

  StatusCode execute(const EventContext&) const override {
    try {
      const auto* input  = m_inputHandle->get();
      auto*       output = m_outputHandle->createAndPut();

      const auto& zpos   = m_zPositions.value();
      const int   nz     = static_cast<int>(zpos.size());
      const long long evtNum  = m_eventCount.fetch_add(1);
      const bool      doPrint = (evtNum % m_debugFreq == 0);

      for (const auto& hit : *input) {
        const std::uint64_t cellID = hit.getCellID();
        const int layer = static_cast<int>(m_decoder->get(cellID, "layer"));

        float newZ;
        if (layer < 0 || layer >= nz) {
          warning() << "[DetectorFlipper] Layer " << layer
                    << " out of range [0," << nz
                    << "); keeping original z." << endmsg;
          newZ = hit.getPosition().z;
        } else {
          newZ = static_cast<float>(zpos[layer]);
        }

        auto nh = output->create();
        nh.setCellID(cellID);
        nh.setEnergy(hit.getEnergy());
        const auto& p = hit.getPosition();
        nh.setPosition({p.x, p.y, newZ});

        if (doPrint) {
          debug() << "  layer=" << layer
                  << "  z: " << p.z << " -> " << newZ << " mm" << endmsg;
        }
      }

      if (doPrint) {
        info() << "DetectorFlipper [evt " << evtNum << "]: "
               << input->size() << " hits remapped." << endmsg;
      }
      return StatusCode::SUCCESS;
    } catch (const std::exception& e) {
      error() << "[DetectorFlipper] Exception in execute(): " << e.what() << endmsg;
      return StatusCode::FAILURE;
    } catch (...) {
      error() << "[DetectorFlipper] Unknown exception in execute()." << endmsg;
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
      error() << "[DetectorFlipper] Exception in finalize(): " << e.what() << endmsg;
      return StatusCode::FAILURE;
    } catch (...) {
      error() << "[DetectorFlipper] Unknown exception in finalize()." << endmsg;
      return StatusCode::FAILURE;
    }
  }

private:
  // Default ZPositions = simulation layer centres extracted from DD4hep geometry
  // (digitized.edm4hep.root, CalorimeterHit.position.z per layer).
  // For real TB data with flipped detector, override with the desired z convention.
  Gaudi::Property<std::string> m_inputName{
      this, "InputCollection", "SiPadHitsDigi",
      "Input SimCalorimeterHit collection"};
  Gaudi::Property<std::string> m_outputName{
      this, "OutputCollection", "SiPadHitsFlipped",
      "Output SimCalorimeterHit collection with remapped z"};
  Gaudi::Property<std::vector<double>> m_zPositions{
      this, "ZPositions",
      {-116.35, -99.75, -83.15, -66.55, -49.95, -33.35, -16.75, -0.15,
        16.45,   33.05,  49.65,  77.25,  93.85, 110.45, 126.98},
      "Target z position [mm] for each layer (indexed 0..NLayers-1)"};
  Gaudi::Property<std::string> m_bitField{
      this, "BitField", "system:8,layer:8,slice:5,x:9,y:9",
      "DD4hep CellID encoding string used to decode the 'layer' field"};
  Gaudi::Property<int> m_nLayers{
      this, "NLayers", 15,
      "Expected number of layers (used to validate ZPositions; 0 = skip)"};
  Gaudi::Property<int> m_debugFreq{
      this, "DebugFrequency", 500, "Print per-event debug info every N events"};

  mutable std::unique_ptr<dd4hep::DDSegmentation::BitFieldCoder> m_decoder;
  mutable std::unique_ptr<k4FWCore::DataHandle<edm4hep::SimCalorimeterHitCollection>> m_inputHandle;
  mutable std::unique_ptr<k4FWCore::DataHandle<edm4hep::SimCalorimeterHitCollection>> m_outputHandle;
  mutable std::atomic<long long> m_eventCount{0};
};

DECLARE_COMPONENT(DetectorFlipper)

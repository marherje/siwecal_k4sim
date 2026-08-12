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
 * The pipeline's job3 steering files read that array from
 * mappings/slab_z_positions.yml, which is the canonical test-beam frame and is
 * shared with the event viewer and k4SiWEcalReco. Never write the array out by
 * hand in a steering file.
 *
 * Where the z comes from
 * ----------------------
 * Leaving ZPositions unset takes the layer z straight from the DD4hep geometry
 * and emits a WARNING saying so. That is a no-op flip -- the hits keep the
 * simulation frame -- which is what you want for a geometry cross-check and not
 * what you want when producing test-beam-frame output, hence the warning.
 *
 * There is deliberately no built-in table. This algorithm used to default to a
 * hardcoded 15-entry array with a 16.6 mm pitch, described in this comment as
 * "the DD4hep z values". It was not: the pitch went 11 mm -> 15 mm in July 2026
 * and the array tracked neither. A wrong z table does not fail, it produces
 * perfectly plausible hits in the wrong place, so the only safe default is the
 * geometry itself.
 *
 * When ZPositions *is* given, its layer spacing is cross-checked against the
 * geometry (|z[i+1]-z[i]|, so a legitimate sign flip or offset does not trip
 * it) and a mismatch is warned about. That is the check that would have caught
 * the stale table.
 *
 * Properties
 * ----------
 *   InputCollection   CalorimeterHit collection to read (default: SiPadHitsDigi)
 *   OutputCollection  CalorimeterHit collection to write (default: SiPadHitsFlipped)
 *   ZPositions        Per-layer z values [mm]. Empty (default) = read DD4hep and
 *                     warn. Otherwise must have exactly NLayers entries.
 *   CompactFile       Compact XML to read the geometry from; may be empty if the
 *                     job already loaded one.
 *   DetectorName      DD4hep detector to take the layers from (default: SiPad).
 *   SpacingToleranceMm  Tolerance of the ZPositions-vs-geometry check (0.5 mm).
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

// Shared with ACTSGeoSvc: one walk of the DD4hep tree, one answer to "where is
// layer N". See SiPadLayerGeometry.h.
#include "SiPadLayerGeometry.h"
#include "DD4hep/Detector.h"

#include <atomic>
#include <cmath>
#include <memory>
#include <sstream>
#include <string>
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

      // ZPositions unset -> take the layer z straight from the DD4hep
      // geometry. There is deliberately no built-in table: this algorithm used
      // to default to a hardcoded array that went stale when the layer pitch
      // changed from 11 mm to 15 mm, and a wrong z table produces perfectly
      // plausible-looking hits. The geometry cannot go stale, because it is
      // what the simulation actually used.
      if (m_zPositions.value().empty()) {
        if (!loadZFromGeometry()) return StatusCode::FAILURE;
      } else {
        m_z = m_zPositions.value();
        info() << "[DetectorFlipper] Using the configured ZPositions ("
               << m_z.size() << " layers)." << endmsg;
        // A supplied table is a second copy of the geometry, so check it still
        // describes the same detector rather than letting it drift unnoticed.
        crossCheckAgainstGeometry();
      }

      const auto& zpos = m_z;
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

  /// Fill m_z from DD4hep and say so loudly. Loud on purpose: it means nobody
  /// declared a target frame, so the hits keep the simulation z and this
  /// algorithm is a no-op — fine for a geometry check, wrong for producing
  /// test-beam-frame output, and the log is the only place that difference
  /// shows up.
  bool loadZFromGeometry() {
    auto& desc = dd4hep::Detector::getInstance();
    try {
      if (!m_compactFile.value().empty()) desc.fromXML(m_compactFile.value());
      m_z = sipad::layerZ(desc, m_detectorName.value());
    } catch (const std::exception& e) {
      error() << "[DetectorFlipper] ZPositions is not set and the layer z could "
                 "not be read from the DD4hep geometry: " << e.what()
              << ". Set CompactFile to the compact XML, or set ZPositions "
                 "explicitly (e.g. from mappings/slab_z_positions.yml)."
              << endmsg;
      return false;
    }
    if (m_z.empty()) {
      error() << "[DetectorFlipper] The DD4hep geometry reports no SiPad layers."
              << endmsg;
      return false;
    }
    std::ostringstream os;
    for (std::size_t i = 0; i < m_z.size(); ++i)
      os << (i ? ", " : "") << m_z[i];
    warning() << "[DetectorFlipper] ZPositions not set: taking the layer z "
                 "directly from the DD4hep geometry"
              << (m_compactFile.value().empty()
                      ? std::string(" already loaded in this job")
                      : " (" + m_compactFile.value() + ")")
              << ". The hits therefore keep the simulation frame and no flip is "
                 "applied. Set ZPositions (e.g. from "
                 "mappings/slab_z_positions.yml) to write a different frame. "
                 "z [mm] = [" << os.str() << "]" << endmsg;
    return true;
  }

  /// Warn when a configured ZPositions no longer matches the layer spacing of
  /// the geometry. Compares |z[i+1]-z[i]|, so it is insensitive to the sign and
  /// offset a frame change legitimately introduces, and catches the thing that
  /// actually went wrong before: a table left at the old pitch. Silent when no
  /// geometry is available to compare against.
  void crossCheckAgainstGeometry() {
    std::vector<double> geo;
    try {
      auto& desc = dd4hep::Detector::getInstance();
      if (!m_compactFile.value().empty()) desc.fromXML(m_compactFile.value());
      geo = sipad::layerZ(desc, m_detectorName.value());
    } catch (const std::exception&) {
      return;  // no geometry loaded in this job; nothing to compare against
    }
    if (geo.size() != m_z.size()) {
      warning() << "[DetectorFlipper] ZPositions has " << m_z.size()
                << " entries but the DD4hep geometry has " << geo.size()
                << " layers." << endmsg;
      return;
    }
    for (std::size_t i = 0; i + 1 < geo.size(); ++i) {
      const double dCfg = std::abs(m_z[i + 1] - m_z[i]);
      const double dGeo = std::abs(geo[i + 1] - geo[i]);
      if (std::abs(dCfg - dGeo) > m_spacingTolerance.value()) {
        warning() << "[DetectorFlipper] ZPositions disagrees with the geometry: "
                     "layers " << i << "->" << (i + 1) << " are " << dCfg
                  << " mm apart in ZPositions but " << dGeo
                  << " mm apart in DD4hep. One of the two is stale; the "
                     "geometry is the one the simulation used." << endmsg;
        return;  // one warning is enough to send someone looking
      }
    }
    info() << "[DetectorFlipper] ZPositions layer spacing matches the DD4hep "
              "geometry." << endmsg;
  }

  StatusCode execute(const EventContext&) const override {
    try {
      const auto* input  = m_inputHandle->get();
      auto*       output = m_outputHandle->createAndPut();

      const auto& zpos   = m_z;
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
      this, "ZPositions", {},
      "Target z position [mm] for each layer (indexed 0..NLayers-1). Empty (the "
      "default) takes the layer z from the DD4hep geometry and warns that it is "
      "doing so, which leaves the hits in the simulation frame. There is no "
      "built-in table on purpose: the one that used to be here went stale when "
      "the layer pitch changed from 11 mm to 15 mm."};
  Gaudi::Property<std::string> m_compactFile{
      this, "CompactFile", "",
      "Compact XML to read the layer z from when ZPositions is empty. May be "
      "left empty if the job already loaded a DD4hep geometry."};
  Gaudi::Property<std::string> m_detectorName{
      this, "DetectorName", "SiPad", "DD4hep detector to take the layers from"};
  Gaudi::Property<double> m_spacingTolerance{
      this, "SpacingToleranceMm", 0.5,
      "How far a configured ZPositions layer spacing may differ from the "
      "geometry before it is flagged as stale [mm]"};
  Gaudi::Property<std::string> m_bitField{
      this, "BitField", "system:8,layer:8,slice:5,x:9,y:9",
      "DD4hep CellID encoding string used to decode the 'layer' field"};
  Gaudi::Property<int> m_nLayers{
      this, "NLayers", 15,
      "Expected number of layers (used to validate ZPositions; 0 = skip)"};
  Gaudi::Property<int> m_debugFreq{
      this, "DebugFrequency", 500, "Print per-event debug info every N events"};

  /// The z table actually in use: either ZPositions or the DD4hep geometry.
  std::vector<double> m_z;

  mutable std::unique_ptr<dd4hep::DDSegmentation::BitFieldCoder> m_decoder;
  mutable std::unique_ptr<k4FWCore::DataHandle<edm4hep::SimCalorimeterHitCollection>> m_inputHandle;
  mutable std::unique_ptr<k4FWCore::DataHandle<edm4hep::SimCalorimeterHitCollection>> m_outputHandle;
  mutable std::atomic<long long> m_eventCount{0};
};

DECLARE_COMPONENT(DetectorFlipper)

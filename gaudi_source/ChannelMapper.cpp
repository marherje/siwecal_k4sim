/*
 * ChannelMapper: translate simulation CellIDs to real TB channel IDs and apply
 * MIP-calibration masking.
 *
 * Motivation
 * ----------
 * The simulation encodes hits with a Cartesian grid CellID
 * ("system:8,layer:8,slice:5,x:9,y:9"), while the real test-beam data uses a
 * per-readout-chip encoding ("system:8,slab:8,chip:16,channel:8,sca:8").
 * Downstream tools (EcalPidTransformer, event viewer, ecal TTree converter)
 * expect the TB encoding so they can cross-check simulation and data directly.
 *
 * This algorithm does two things per hit:
 *   1. CellID remapping — finds the nearest physical pad in the FEV10/FEV11 pad
 *      map (nearest-neighbour in (x, y)) and rewrites the CellID to
 *      system:8,slab:8,chip:16,channel:8,sca:8  (sca = 0 for simulation).
 *   2. MIP masking — checks the remapped (slab, chip, channel) against the MIP
 *      calibration file.  Channels where mpv == 0, mpv is NaN, or mpv > MaxMIPValue
 *      are flagged.  Hits outside the pad-map range (no nearest neighbour within
 *      PositionTolerance mm) are also flagged.
 *
 * Outputs (written in parallel, index-aligned):
 *   OutputCollection   — SimCalorimeterHitCollection with new CellIDs; positions
 *                        and energies are copied unchanged.
 *   OutputMaskedFlags  — podio::UserDataCollection<int32_t>, one entry per hit:
 *                        0 = calibrated channel, 1 = masked / unmapped.
 *
 * The masking flag flows into analysis/sim_to_ecal_tree.py which writes it as
 * the hit_ismasked branch in the ecal TTree consumed by EcalToEDM4hep /
 * EcalPidTransformer.
 *
 * Pad map format  (fev10_rotate_chip_channel_x_y_mapping.txt)
 * -----------------------------------------------------------
 *   chip  x0  y0  channel  x  y
 *   0  67.05  67.05  16  86.3  86.3
 *   ...
 * (first non-comment line is a text header, silently skipped; 1024 data rows,
 * 16 chips × 64 channels; x/y in mm)
 *
 * MIP calibration format  (dummy_mip_map_15_highgain.txt or real files)
 * -----------------------------------------------------------------------
 *   # comment
 *   layer  chip  channel  mpv  empv  widthmpv  chi2ndf  nentries
 *   0  0  0  20.0  0.0  2.0  1.0  1000
 *
 * Properties
 * ----------
 *   InputCollection      SimCalorimeterHit collection to read (default: SiPadHitsFlipped)
 *   OutputCollection     Remapped hit collection name  (default: SiPadHitsMapped)
 *   OutputMaskedFlags    Masking flag collection name   (default: SiPadHitsMasked)
 *   PadMapFile           Path to default FEV10 pad map
 *   PadMapFileSlab12     Path to FEV11 pad map for slab 12 (empty = use default)
 *   MIPCalibFile         Path to merged MIP calibration file
 *   MaxMIPValue          Channels with mpv > this are masked (default: 100.0)
 *   PositionTolerance    Max (x,y) distance [mm] to accept a pad match (default: 4.0)
 *   BitFieldIn           Input CellID encoding string
 *   BitFieldOut          Output CellID encoding string
 *   DebugFrequency       Print per-event debug info every N events
 */
#include "k4FWCore/Transformer.h"

#include "edm4hep/SimCalorimeterHitCollection.h"
#include "podio/UserDataCollection.h"

#include "DDSegmentation/BitFieldCoder.h"

#include "Gaudi/Property.h"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <memory>
#include <set>
#include <sstream>
#include <string>
#include <tuple>
#include <vector>

// --------------------------------------------------------------------------- //
using OutputType = std::tuple<edm4hep::SimCalorimeterHitCollection,
                              podio::UserDataCollection<std::int32_t>>;

struct ChannelMapper final
    : k4FWCore::MultiTransformer<OutputType(const edm4hep::SimCalorimeterHitCollection&)> {

  ChannelMapper(const std::string& name, ISvcLocator* svc)
      : MultiTransformer(
            name, svc,
            {KeyValues("InputCollection", {"SiPadHitsFlipped"})},
            {KeyValues("OutputCollection", {"SiPadHitsMapped"}),
             KeyValues("OutputMaskedFlags", {"SiPadHitsMasked"})}) {}

  // ----------------------------------------------------------------- //
  // initialize: load pad map(s) and MIP calibration
  // ----------------------------------------------------------------- //
  StatusCode initialize() override {
    StatusCode sc = MultiTransformer::initialize();
    if (sc.isFailure()) return sc;

    // --- input bitfield decoder ---
    if (m_bitFieldIn.value().empty()) {
      error() << "[ChannelMapper] BitFieldIn must be set." << endmsg;
      return StatusCode::FAILURE;
    }
    m_decoderIn = std::make_unique<dd4hep::DDSegmentation::BitFieldCoder>(m_bitFieldIn.value());

    // --- output bitfield coder ---
    if (m_bitFieldOut.value().empty()) {
      error() << "[ChannelMapper] BitFieldOut must be set." << endmsg;
      return StatusCode::FAILURE;
    }
    m_decoderOut = std::make_unique<dd4hep::DDSegmentation::BitFieldCoder>(m_bitFieldOut.value());

    // --- load pad maps ---
    if (!loadPadMap(m_padMapFile.value(), m_padDefault)) return StatusCode::FAILURE;
    info() << "[ChannelMapper] Default pad map: " << m_padDefault.size()
           << " entries from " << m_padMapFile.value() << endmsg;

    if (!m_padMapFileSlab12.value().empty()) {
      if (!loadPadMap(m_padMapFileSlab12.value(), m_padSlab12)) return StatusCode::FAILURE;
      info() << "[ChannelMapper] Slab-12 pad map: " << m_padSlab12.size()
             << " entries from " << m_padMapFileSlab12.value() << endmsg;
    }

    // --- load MIP calibration ---
    if (!loadMIPCalib(m_mipCalibFile.value())) return StatusCode::FAILURE;
    info() << "[ChannelMapper] MIP calibration: " << m_masked.size()
           << " masked (slab,chip,chan) triples from " << m_mipCalibFile.value() << endmsg;

    return StatusCode::SUCCESS;
  }

  // ----------------------------------------------------------------- //
  // finalize
  // ----------------------------------------------------------------- //
  StatusCode finalize() override {
    m_decoderIn.reset();
    m_decoderOut.reset();
    return MultiTransformer::finalize();
  }

  // ----------------------------------------------------------------- //
  // operator(): one event → (remapped hits, masking flags)
  // ----------------------------------------------------------------- //
  OutputType operator()(const edm4hep::SimCalorimeterHitCollection& input) const override {
    edm4hep::SimCalorimeterHitCollection outColl;
    podio::UserDataCollection<std::int32_t> maskedFlags;

    const long long evtNum = m_eventCount.fetch_add(1);
    const bool      doPrint = (evtNum % m_debugFreq.value() == 0);
    int nMasked = 0;

    const auto& pads = m_padDefault;  // alias for readability

    for (const auto& hit : input) {
      const std::uint64_t cellID = hit.getCellID();
      const int slab  = static_cast<int>(m_decoderIn->get(cellID, "layer"));
      const auto& pos = hit.getPosition();
      const float hx  = pos.x, hy = pos.y;

      // --- nearest-neighbour lookup in the pad map ---
      const auto& padMap = (!m_padSlab12.empty() && slab == 12) ? m_padSlab12 : pads;

      int   bestChip = 0, bestChan = 0;
      float bestDist2 = 1e12f;
      for (const auto& p : padMap) {
        const float dx = p.x - hx, dy = p.y - hy;
        const float d2 = dx*dx + dy*dy;
        if (d2 < bestDist2) { bestDist2 = d2; bestChip = p.chip; bestChan = p.channel; }
      }

      const float tol    = static_cast<float>(m_posTolerance.value());
      const bool unmapped = std::sqrt(bestDist2) > tol;

      // --- build output CellID (TB encoding) ---
      std::uint64_t newCellID = 0;
      m_decoderOut->set(newCellID, "system", m_decoderIn->get(cellID, "system"));
      m_decoderOut->set(newCellID, "slab",    slab);
      if (!unmapped) {
        m_decoderOut->set(newCellID, "chip",    bestChip);
        m_decoderOut->set(newCellID, "channel", bestChan);
      }
      // sca = 0 (no switched-capacitor array concept in simulation)

      // --- masking check ---
      const bool isMasked = unmapped ||
                            m_masked.count(std::make_tuple(slab, bestChip, bestChan)) > 0;

      // --- write output ---
      auto nh = outColl.create();
      nh.setCellID(newCellID);
      nh.setEnergy(hit.getEnergy());
      nh.setPosition(pos);

      maskedFlags.push_back(isMasked ? 1 : 0);

      if (isMasked) ++nMasked;
    }

    if (doPrint) {
      info() << "ChannelMapper [evt " << evtNum << "]: "
             << input.size() << " hits, " << nMasked << " masked." << endmsg;
    }
    return {std::move(outColl), std::move(maskedFlags)};
  }

private:
  // ----------------------------------------------------------------- //
  // Helpers
  // ----------------------------------------------------------------- //
  struct PadEntry { float x, y; int chip, channel; };

  bool loadPadMap(const std::string& path, std::vector<PadEntry>& dest) const {
    std::ifstream fin(path);
    if (!fin) {
      error() << "[ChannelMapper] Cannot open pad map file: " << path << endmsg;
      return false;
    }
    std::string line;
    while (std::getline(fin, line)) {
      if (line.empty() || line[0] == '#') continue;
      std::istringstream ss(line);
      int chip, channel;
      float x0, y0, x, y;
      // Format: chip x0 y0 channel x y
      // Skip the text header ("chip x0 y0 ...") by checking first token is numeric
      ss >> chip;
      if (ss.fail()) continue;  // non-numeric first token → header row
      ss >> x0 >> y0 >> channel >> x >> y;
      if (ss.fail()) continue;
      dest.push_back({x, y, chip, channel});
    }
    if (dest.empty()) {
      error() << "[ChannelMapper] Pad map loaded 0 entries from " << path << endmsg;
      return false;
    }
    return true;
  }

  bool loadMIPCalib(const std::string& path) const {
    std::ifstream fin(path);
    if (!fin) {
      error() << "[ChannelMapper] Cannot open MIP calibration file: " << path << endmsg;
      return false;
    }
    const float maxMIP = static_cast<float>(m_maxMIP.value());
    std::string line;
    int nRead = 0;
    while (std::getline(fin, line)) {
      if (line.empty() || line[0] == '#') continue;
      std::istringstream ss(line);
      int layer, chip, chan;
      float mpv;
      ss >> layer >> chip >> chan >> mpv;
      if (ss.fail()) continue;
      ++nRead;
      // Mask if mpv == 0, NaN, or above the maximum credible MIP value
      if (mpv == 0.0f || std::isnan(mpv) || mpv > maxMIP) {
        m_masked.insert(std::make_tuple(layer, chip, chan));
      }
    }
    if (nRead == 0) {
      error() << "[ChannelMapper] MIP calibration file has no data rows: " << path << endmsg;
      return false;
    }
    return true;
  }

  // ----------------------------------------------------------------- //
  // Properties
  // ----------------------------------------------------------------- //
  Gaudi::Property<std::string> m_padMapFile{
      this, "PadMapFile",
      "masking_info/geometry/fev10_rotate_chip_channel_x_y_mapping.txt",
      "Path to FEV10 pad map (chip/channel → x,y) for the default slab"};
  Gaudi::Property<std::string> m_padMapFileSlab12{
      this, "PadMapFileSlab12", "",
      "Path to FEV11 pad map for slab 12 override (empty = use default map)"};
  Gaudi::Property<std::string> m_mipCalibFile{
      this, "MIPCalibFile",
      "masking_info/calibration/dummy_mip_map_15_highgain.txt",
      "Path to merged MIP calibration file (layer chip channel mpv ...)"};
  Gaudi::Property<double> m_maxMIP{
      this, "MaxMIPValue", 100.0,
      "Channels with mpv > MaxMIPValue are masked (typical: 100 MIP)"};
  Gaudi::Property<double> m_posTolerance{
      this, "PositionTolerance", 4.0,
      "Maximum Euclidean distance [mm] from a hit to the nearest pad centre "
      "before the hit is considered unmapped (and masked).  The default 4.0 mm "
      "covers the half-gap at FEV10 chip boundaries (~3.3 mm)."};
  Gaudi::Property<std::string> m_bitFieldIn{
      this, "BitFieldIn", "system:8,layer:8,slice:5,x:9,y:9",
      "DD4hep CellID encoding of the input collection"};
  Gaudi::Property<std::string> m_bitFieldOut{
      this, "BitFieldOut", "system:8,slab:8,chip:16,channel:8,sca:8",
      "DD4hep CellID encoding of the output collection"};
  Gaudi::Property<int> m_debugFreq{
      this, "DebugFrequency", 500, "Print per-event summary every N events"};

  // ----------------------------------------------------------------- //
  // State initialised in initialize()
  // ----------------------------------------------------------------- //
  std::vector<PadEntry> m_padDefault;
  std::vector<PadEntry> m_padSlab12;
  mutable std::set<std::tuple<int,int,int>> m_masked;

  std::unique_ptr<dd4hep::DDSegmentation::BitFieldCoder> m_decoderIn;
  std::unique_ptr<dd4hep::DDSegmentation::BitFieldCoder> m_decoderOut;

  mutable std::atomic<long long> m_eventCount{0};
};

DECLARE_COMPONENT(ChannelMapper)

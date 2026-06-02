#include "Gaudi/Algorithm.h"
#include "GaudiKernel/MsgStream.h"
#include "k4FWCore/DataHandle.h"
#include "edm4hep/SimCalorimeterHitCollection.h"
#include "edm4hep/TrackerHit3DCollection.h"
#include "edm4hep/TrackCollection.h"
#include "DD4hep/BitFieldCoder.h"
#include "podio/ROOTReader.h"
#include "podio/Frame.h"
#include "ROOT/RNTupleModel.hxx"
#include "ROOT/RNTupleWriter.hxx"
#include "TFile.h"
#include "ValidationUtils.h"

#include <algorithm>
#include <atomic>
#include <cstdint>
#include <map>
#include <memory>
#include <optional>
#include <string>
#include <vector>

class EDM4HEP2RNTuple : public Gaudi::Algorithm {
public:
  EDM4HEP2RNTuple(const std::string& name, ISvcLocator* svcLoc)
      : Gaudi::Algorithm(name, svcLoc) {}

  StatusCode initialize() override {
    try {
      StatusCode sc = Gaudi::Algorithm::initialize();
      if (sc.isFailure()) return sc;

      const auto& colls       = m_collections.value();
      const auto& bfs         = m_bitFields.value();
      const auto& params      = m_sourceIDParams.value();
      const auto& detIDs      = m_detectorIDs.value();
      const auto& ntupleNames = m_ntupleNames.value();

      // 1. Validate parallel vectors
      if (colls.size() != bfs.size()         ||
          colls.size() != params.size()      ||
          colls.size() != detIDs.size()      ||
          colls.size() != ntupleNames.size()) {
        error() << "[EDM4HEP2RNTuple] Collections, BitFields, SourceIDParams,"
                << " DetectorIDs and NTupleNames must all have the same size."
                << " Got: "
                << colls.size()       << " " << bfs.size()   << " "
                << params.size()      << " " << detIDs.size() << " "
                << ntupleNames.size()
                << endmsg;
        return StatusCode::FAILURE;
      }

      // 2. Build one DataHandle per collection
      for (const auto& cname : colls) {
        m_inputHandles.push_back(
            std::make_unique<k4FWCore::DataHandle<edm4hep::SimCalorimeterHitCollection>>(
                cname, Gaudi::DataHandle::Reader, this));
      }

      // 3. Build one BitFieldCoder per collection
      try {
        for (const auto& bf : bfs) {
          m_decoders.push_back(
              std::make_unique<dd4hep::DDSegmentation::BitFieldCoder>(bf));
        }
      } catch (const std::exception& e) {
        error() << "[EDM4HEP2RNTuple] Invalid BitField string: " << e.what() << endmsg;
        return StatusCode::FAILURE;
      }

      // 4. Inspect BitField field names per detector
      for (size_t i = 0; i < m_decoders.size(); ++i) {
        const auto& decoder = *m_decoders[i];
        std::vector<std::string> names;
        for (size_t f = 0; f < decoder.size(); ++f) {
          names.push_back(decoder[f].name());
        }
        m_detectorFieldNames.push_back(names);

        info() << "[EDM4HEP2RNTuple] detector_id=" << detIDs[i]
               << " ntuple='" << ntupleNames[i] << "' BitField fields: ";
        for (const auto& n : names) info() << n << " ";
        info() << endmsg;
      }

      // 5. Open podio::ROOTReader for frame parameter access
      m_reader = std::make_unique<podio::ROOTReader>();
      m_reader->openFile(m_inputFile.value());

      if (m_reader->getEntries("events") == 0) {
        error() << "[EDM4HEP2RNTuple] Input file has 0 events." << endmsg;
        return StatusCode::FAILURE;
      }

      // Validate first frame
      {
        auto firstFrame = podio::Frame(m_reader->readEntry("events", 0));
        for (const auto& cname : colls)
          validateCollection(firstFrame, cname, "EDM4HEP2RNTuple", msgStream());
        for (const auto& pname : params) {
          if (!firstFrame.getParameter<std::vector<int>>(pname).has_value()) {
            warning() << "[EDM4HEP2RNTuple] Frame parameter '" << pname
                      << "' not found in first window. source_id will be 0." << endmsg;
          }
        }
      }

      // 6. Open output TFile once; all per-detector writers share it
      m_outputTFile.reset(
          TFile::Open(m_outputFile.value().c_str(), "RECREATE"));
      if (!m_outputTFile || m_outputTFile->IsZombie()) {
        error() << "[EDM4HEP2RNTuple] Could not open output file: "
                << m_outputFile.value() << endmsg;
        return StatusCode::FAILURE;
      }

      // 7. Build one RNTuple model and writer per detector
      for (size_t i = 0; i < colls.size(); ++i) {
        auto model = ROOT::RNTupleModel::Create();
        DetectorFields df;

        // Fixed fields (no detector_id: each RNTuple is already detector-specific)
        auto fWindowID     = model->MakeField<int>("window_id");
        auto fCellID       = model->MakeField<uint64_t>("cellID");
        auto fX            = model->MakeField<float>("x");
        auto fY            = model->MakeField<float>("y");
        auto fZ            = model->MakeField<float>("z");
        auto fE            = model->MakeField<float>("E");
        auto fT            = model->MakeField<float>("t");
        auto fSourceID     = model->MakeField<int>("source_id");
        auto fTWindowStart = model->MakeField<float>("t_window_start");

        df.fWindowID     = fWindowID.get();
        df.fCellID       = fCellID.get();
        df.fX            = fX.get();
        df.fY            = fY.get();
        df.fZ            = fZ.get();
        df.fE            = fE.get();
        df.fT            = fT.get();
        df.fSourceID     = fSourceID.get();
        df.fTWindowStart = fTWindowStart.get();

        // Dynamic BitField columns for THIS detector only — no cross-detector padding
        info() << "[EDM4HEP2RNTuple] NTuple '" << ntupleNames[i]
               << "' dynamic BitField columns: ";
        for (const auto& fname : m_detectorFieldNames[i]) {
          const std::string colName = "bf_" + fname;
          auto fptr = model->MakeField<int>(colName);
          df.fieldPtrs[fname] = fptr.get();
          df.fieldSharedPtrs.push_back(std::move(fptr));
          info() << colName << " ";
        }
        info() << endmsg;

        // File was opened empty with "RECREATE"; Append works for all detectors.
        m_writers.push_back(
            ROOT::RNTupleWriter::Append(
                std::move(model), ntupleNames[i], *m_outputTFile));

        m_detFields.push_back(std::move(df));

        info() << "[EDM4HEP2RNTuple]   detector_id=" << detIDs[i]
               << "  collection=" << colls[i]
               << "  ntuple=" << ntupleNames[i]
               << "  bitfield=" << bfs[i]
               << endmsg;
      }

      // ---- Measurement (TrackerHit3D) RNTuple setup -------------------------
      const auto& mColls   = m_measCollections.value();
      const auto& mNames   = m_measNtupleNames.value();
      const auto& mBFields = m_measBitFields.value();

      if (!mColls.empty()) {
        if (mColls.size() != mNames.size() || mColls.size() != mBFields.size()) {
          error() << "[EDM4HEP2RNTuple] MeasCollections, MeasNtupleNames and "
                  << "MeasBitFields must have the same size." << endmsg;
          return StatusCode::FAILURE;
        }

        for (std::size_t mi = 0; mi < mColls.size(); ++mi) {
          m_measHandles.push_back(
              std::make_unique<k4FWCore::DataHandle<edm4hep::TrackerHit3DCollection>>(
                  mColls[mi], Gaudi::DataHandle::Reader, this));

          m_measDecoders.push_back(
              std::make_unique<dd4hep::DDSegmentation::BitFieldCoder>(mBFields[mi]));

          auto model = ROOT::RNTupleModel::Create();

          auto fWindowID = model->MakeField<int>("window_id");
          auto fX        = model->MakeField<float>("x");
          auto fY        = model->MakeField<float>("y");
          auto fZ        = model->MakeField<float>("z");
          auto fPlane    = model->MakeField<int>("plane");
          auto fLayer    = model->MakeField<int>("layer");
          auto fColumn   = model->MakeField<int>("column");
          auto fRow      = model->MakeField<int>("row");
          auto fEdep     = model->MakeField<float>("edep");
          auto fSourceID = model->MakeField<int>("source_id");

          MeasFields mf;
          mf.fWindowID = fWindowID.get();
          mf.fX        = fX.get();
          mf.fY        = fY.get();
          mf.fZ        = fZ.get();
          mf.fPlane    = fPlane.get();
          mf.fLayer    = fLayer.get();
          mf.fColumn   = fColumn.get();
          mf.fRow      = fRow.get();
          mf.fEdep     = fEdep.get();
          mf.fSourceID = fSourceID.get();

          mf.intPtrs.push_back(std::move(fWindowID));
          mf.intPtrs.push_back(std::move(fPlane));
          mf.intPtrs.push_back(std::move(fLayer));
          mf.intPtrs.push_back(std::move(fColumn));
          mf.intPtrs.push_back(std::move(fRow));
          mf.intPtrs.push_back(std::move(fSourceID));
          mf.floatPtrs.push_back(std::move(fX));
          mf.floatPtrs.push_back(std::move(fY));
          mf.floatPtrs.push_back(std::move(fZ));
          mf.floatPtrs.push_back(std::move(fEdep));

          m_measFields.push_back(std::move(mf));

          m_measWriters.push_back(ROOT::RNTupleWriter::Append(
              std::move(model), mNames[mi], *m_outputTFile));

          info() << "[EDM4HEP2RNTuple] Measurement RNTuple '" << mNames[mi]
                 << "' from collection '" << mColls[mi] << "'" << endmsg;
        }
      }

      // ---- Track RNTuple setup (optional) ----------------------------------
      if (!m_trackFile.value().empty()) {
        try {
          // Open track file with direct podio reader (bypass Gaudi DataHandle)
          m_trackReader = std::make_unique<podio::ROOTReader>();
          m_trackReader->openFile(m_trackFile.value());
          m_nTrackEvents = m_trackReader->getEntries("events");

          // Validate event count matches hit file
          const std::size_t nHitEvents =
              static_cast<std::size_t>(m_reader->getEntries("events"));
          if (m_nTrackEvents != nHitEvents) {
            error() << "[EDM4HEP2RNTuple] Track file '"
                    << m_trackFile.value() << "' has "
                    << m_nTrackEvents << " event(s) but hit file '"
                    << m_inputFile.value() << "' has " << nHitEvents
                    << " event(s). Both must have the same number of events."
                    << endmsg;
            return StatusCode::FAILURE;
          }

          // Build minimal track RNTuple: one row per track
          auto model = ROOT::RNTupleModel::Create();

          auto fWindowID = model->MakeField<int>("window_id");
          auto fTrackID  = model->MakeField<int>("track_id");
          auto fChi2     = model->MakeField<float>("chi2");
          auto fNdf      = model->MakeField<int>("ndf");
          auto fSeedX    = model->MakeField<float>("seed_x");  // DD4hep transverse X [mm]
          auto fSeedY    = model->MakeField<float>("seed_y");  // DD4hep transverse Y [mm]

          m_fTrackWindowID = fWindowID.get();
          m_fTrackID       = fTrackID.get();
          m_fTrackChi2     = fChi2.get();
          m_fTrackNdf      = fNdf.get();
          m_fTrackSeedX    = fSeedX.get();
          m_fTrackSeedY    = fSeedY.get();

          m_trackIntPtrs.push_back(std::move(fWindowID));
          m_trackIntPtrs.push_back(std::move(fTrackID));
          m_trackIntPtrs.push_back(std::move(fNdf));
          m_trackFloatPtrs.push_back(std::move(fChi2));
          m_trackFloatPtrs.push_back(std::move(fSeedX));
          m_trackFloatPtrs.push_back(std::move(fSeedY));

          m_trackWriter = ROOT::RNTupleWriter::Append(
            std::move(model), m_trackCollectionName.value(), *m_outputTFile);

          // Per-surface track states RNTuple
          {
            auto tsModel  = ROOT::RNTupleModel::Create();
            auto fWin     = tsModel->MakeField<int>("window_id");
            auto fTrkId   = tsModel->MakeField<int>("track_id");
            auto fStIdx   = tsModel->MakeField<int>("state_idx");
            auto fX       = tsModel->MakeField<float>("x");
            auto fY       = tsModel->MakeField<float>("y");
            auto fZ       = tsModel->MakeField<float>("z");

            auto fLoc0  = tsModel->MakeField<float>("loc0");
            auto fTilt  = tsModel->MakeField<float>("tilt");
            auto fChi2ts = tsModel->MakeField<float>("chi2");
            auto fNmeas  = tsModel->MakeField<int>("nmeas");

            m_fTsWindowID = fWin.get();
            m_fTsTrackID  = fTrkId.get();
            m_fTsStateIdx = fStIdx.get();
            m_fTsX        = fX.get();
            m_fTsY        = fY.get();
            m_fTsZ        = fZ.get();
            m_fTsLoc0     = fLoc0.get();
            m_fTsTilt     = fTilt.get();
            m_fTsChi2     = fChi2ts.get();
            m_fTsNmeas    = fNmeas.get();

            m_tsIntPtrs.push_back(std::move(fWin));
            m_tsIntPtrs.push_back(std::move(fTrkId));
            m_tsIntPtrs.push_back(std::move(fStIdx));
            m_tsIntPtrs.push_back(std::move(fNmeas));
            m_tsFloatPtrs.push_back(std::move(fX));
            m_tsFloatPtrs.push_back(std::move(fY));
            m_tsFloatPtrs.push_back(std::move(fZ));
            m_tsFloatPtrs.push_back(std::move(fLoc0));
            m_tsFloatPtrs.push_back(std::move(fTilt));
            m_tsFloatPtrs.push_back(std::move(fChi2ts));

            m_trackStatesWriter = ROOT::RNTupleWriter::Append(
                std::move(tsModel), "ACTSTrackStates", *m_outputTFile);
          }

          info() << "[EDM4HEP2RNTuple] Track export enabled: '"
                 << m_trackFile.value() << "' ("
                 << m_nTrackEvents << " events)" << endmsg;

        } catch (const std::exception& e) {
          error() << "[EDM4HEP2RNTuple] Cannot open track file '"
                  << m_trackFile.value() << "': " << e.what() << endmsg;
          return StatusCode::FAILURE;
        }
      }

      info() << "[EDM4HEP2RNTuple] Initialized."
             << "  Output: " << m_outputFile.value()
             << "  Detectors: " << colls.size() << endmsg;

      return sc;
    } catch (const std::exception& e) {
      error() << "[EDM4HEP2RNTuple] Exception in initialize(): " << e.what() << endmsg;
      return StatusCode::FAILURE;
    } catch (...) {
      error() << "[EDM4HEP2RNTuple] Unknown exception in initialize()." << endmsg;
      return StatusCode::FAILURE;
    }
  }

  // Called once per time window by Gaudi (EvtMax=-1, IOSvc reads the file).
  StatusCode execute(const EventContext&) const override {
    try {
      const long long windowID = m_windowCount.fetch_add(1);
      const auto& colls  = m_collections.value();
      const auto& params = m_sourceIDParams.value();

      // Read the frame once per window to get all parameters
      std::optional<podio::Frame> maybeFrame;
      float t_window_start = 0.f;
      try {
        auto entry = m_reader->readEntry("events", static_cast<size_t>(windowID));
        maybeFrame.emplace(std::move(entry));
        auto ts = maybeFrame->getParameter<float>("t_window_start");
        if (ts.has_value()) t_window_start = ts.value();
      } catch (...) {
        warning() << "[EDM4HEP2RNTuple] Could not read frame parameters for"
                  << " window " << windowID << ". source_id and t_window_start will be 0."
                  << endmsg;
      }

      for (size_t i = 0; i < colls.size(); ++i) {
        const auto* hits = m_inputHandles[i]->get();
        auto& df = m_detFields[i];

        std::vector<int> sourceIDs;
        if (maybeFrame) {
          sourceIDs = maybeFrame->getParameter<std::vector<int>>(params[i])
                          .value_or(std::vector<int>{});
        }

        for (int hitIdx = 0; hitIdx < static_cast<int>(hits->size()); ++hitIdx) {
          const auto& hit       = (*hits)[hitIdx];
          const uint64_t cellID = hit.getCellID();
          const auto& pos       = hit.getPosition();

          float t = 0.f;
          if (hit.contributions_size() > 0)
            t = hit.getContributions(0).getTime();

          const int source_id = (hitIdx < static_cast<int>(sourceIDs.size()))
                                ? sourceIDs[hitIdx] : 0;

          // Decode BitField fields for this detector only (no zero-padding needed)
          for (auto& [fname, ptr] : df.fieldPtrs) {
            *ptr = static_cast<int>(m_decoders[i]->get(cellID, fname));
          }

          *df.fWindowID     = static_cast<int>(windowID);
          *df.fCellID       = cellID;
          *df.fX            = pos.x;
          *df.fY            = pos.y;
          *df.fZ            = pos.z;
          *df.fE            = hit.getEnergy();
          *df.fT            = t;
          *df.fSourceID     = source_id;
          *df.fTWindowStart = t_window_start;

          m_writers[i]->Fill();
        }
      }

      if (windowID % 100 == 0) {
        info() << "[EDM4HEP2RNTuple] Processed window " << windowID
               << "  t_window_start=" << t_window_start << " ns" << endmsg;
        for (size_t i = 0; i < colls.size(); ++i) {
          info() << "[EDM4HEP2RNTuple]   collection=" << colls[i]
                 << "  hits=" << m_inputHandles[i]->get()->size() << endmsg;
        }
      }

      // ---- Export TrackerHit3D measurements ----------------------------------
      for (std::size_t mi = 0; mi < m_measHandles.size(); ++mi) {
        try {
          const auto* hits = m_measHandles[mi]->get();
          if (!hits) continue;

          auto& mf      = m_measFields[mi];
          auto& decoder = *m_measDecoders[mi];
          const int sourceID = 0;

          for (std::size_t hi = 0; hi < hits->size(); ++hi) {
            const auto& hit = (*hits)[hi];
            const auto& pos = hit.getPosition();
            const uint64_t cid = hit.getCellID();

            int plane  = 0, layer = 0, column = 0, row = 0;
            try {
              for (std::size_t fi = 0; fi < decoder.size(); ++fi) {
                const auto& fname = decoder[fi].name();
                if      (fname == "plane")  plane  = static_cast<int>(decoder.get(cid, "plane"));
                else if (fname == "layer")  layer  = static_cast<int>(decoder.get(cid, "layer"));
                else if (fname == "column") column = static_cast<int>(decoder.get(cid, "column"));
                else if (fname == "row")    row    = static_cast<int>(decoder.get(cid, "row"));
              }
            } catch (...) {}

            *mf.fWindowID = static_cast<int>(windowID);
            *mf.fX        = static_cast<float>(pos.x);
            *mf.fY        = static_cast<float>(pos.y);
            *mf.fZ        = static_cast<float>(pos.z);
            *mf.fPlane    = plane;
            *mf.fLayer    = layer;
            *mf.fColumn   = column;
            *mf.fRow      = row;
            *mf.fEdep     = hit.getEDep();
            *mf.fSourceID = sourceID;

            m_measWriters[mi]->Fill();
          }
        } catch (const std::exception& e) {
          warning() << "[EDM4HEP2RNTuple] Could not export measurements["
                    << mi << "] for window " << windowID
                    << ": " << e.what() << " — skipping." << endmsg;
        }
      }

      // ---- Export tracks if enabled ----------------------------------------
      if (m_trackWriter && m_trackReader &&
          static_cast<std::size_t>(windowID) < m_nTrackEvents) {
        try {
          auto entry  = m_trackReader->readEntry("events",
                            static_cast<std::size_t>(windowID));
          auto tframe = podio::Frame(std::move(entry));
          const auto& tracks = tframe.get<edm4hep::TrackCollection>(m_trackCollectionName.value());

          for (int iTrack = 0; iTrack < static_cast<int>(tracks.size()); ++iTrack) {
            const auto& track  = tracks[iTrack];
            const auto& states = track.getTrackStates();

            // Seed position stored in first TrackState (location=AtIP)
            // D0 = seed transverse X [mm], Z0 = seed transverse Y [mm]
            float seedX = 0.0f, seedY = 0.0f;
            for (const auto& ts : states) {
              if (ts.location == edm4hep::TrackState::AtIP) {
                seedX = ts.D0;
                seedY = ts.Z0;
                break;
              }
            }

            *m_fTrackWindowID = static_cast<int>(windowID);
            *m_fTrackID       = iTrack;
            *m_fTrackChi2     = track.getChi2();
            *m_fTrackNdf      = track.getNdf();
            *m_fTrackSeedX    = seedX;
            *m_fTrackSeedY    = seedY;
            m_trackWriter->Fill();

            // Per-surface states (AtOther, tip→seed order; reverse for beam order)
            if (m_trackStatesWriter) {
              std::vector<const edm4hep::TrackState*> surfStates;
              for (const auto& ts : track.getTrackStates()) {
                if (ts.location == edm4hep::TrackState::AtOther)
                  surfStates.push_back(&ts);
              }
              std::reverse(surfStates.begin(), surfStates.end());
              const int nMeasStates = static_cast<int>(surfStates.size());
              for (int si = 0; si < nMeasStates; ++si) {
                const auto& ts  = *surfStates[si];
                *m_fTsWindowID  = static_cast<int>(windowID);
                *m_fTsTrackID   = iTrack;
                *m_fTsStateIdx  = si;
                *m_fTsX         = ts.referencePoint.x;
                *m_fTsY         = ts.referencePoint.y;
                *m_fTsZ         = ts.referencePoint.z;
                *m_fTsLoc0      = ts.D0;     // raw ACTS eBoundLoc0 (strip measurement)
                *m_fTsTilt      = ts.omega;  // surface stereo tilt: ±sin(α) or 0
                *m_fTsChi2      = ts.Z0;     // per-state innovation chi² (carried in Z0)
                *m_fTsNmeas     = nMeasStates;
                m_trackStatesWriter->Fill();
              }
            }
          }
        } catch (const std::exception& e) {
          warning() << "[EDM4HEP2RNTuple] Could not read tracks for window "
                    << windowID << ": " << e.what()
                    << " — skipping." << endmsg;
        }
      }

      return StatusCode::SUCCESS;
    } catch (const std::exception& e) {
      error() << "[EDM4HEP2RNTuple] Exception in execute(): " << e.what() << endmsg;
      return StatusCode::FAILURE;
    } catch (...) {
      error() << "[EDM4HEP2RNTuple] Unknown exception in execute()." << endmsg;
      return StatusCode::FAILURE;
    }
  }

  StatusCode finalize() override {
    try {
      for (auto& w : m_measWriters) if (w) w.reset();
      m_measHandles.clear();
      m_measFields.clear();
      m_measDecoders.clear();

      if (m_trackStatesWriter) { m_trackStatesWriter.reset(); }
      m_tsIntPtrs.clear();
      m_tsFloatPtrs.clear();

      if (m_trackWriter) { m_trackWriter.reset(); }
      m_trackReader.reset();
      m_trackIntPtrs.clear();
      m_trackFloatPtrs.clear();

      for (auto& w : m_writers) w.reset();  // flush all writers before closing file
      m_writers.clear();
      m_detFields.clear();
      m_outputTFile.reset();  // closes the TFile after all writers are flushed
      m_reader.reset();
      for (auto& h : m_inputHandles) h.reset();

      info() << "[EDM4HEP2RNTuple] Done. Total windows processed: "
             << m_windowCount.load() << endmsg;
      info() << "[EDM4HEP2RNTuple] RNTuples written to: " << m_outputFile.value() << endmsg;

      return Gaudi::Algorithm::finalize();
    } catch (const std::exception& e) {
      error() << "[EDM4HEP2RNTuple] Exception in finalize(): " << e.what() << endmsg;
      return StatusCode::FAILURE;
    } catch (...) {
      error() << "[EDM4HEP2RNTuple] Unknown exception in finalize()." << endmsg;
      return StatusCode::FAILURE;
    }
  }

private:
  // Python-configurable properties — parallel vectors, one entry per detector.
  Gaudi::Property<std::vector<std::string>> m_collections{
      this, "Collections",
      {"SiTargetHitsWindowed", "SiPadHitsWindowed"},
      "Input SimCalorimeterHit collection names, one per detector"};
  Gaudi::Property<std::vector<std::string>> m_bitFields{
      this, "BitFields",
      {"system:8,layer:8,slice:4,plane:1,strip:14",
       "system:8,layer:8,slice:4,x:9,y:9"},
      "DD4hep BitField encoding string for each collection"};
  Gaudi::Property<std::vector<std::string>> m_sourceIDParams{
      this, "SourceIDParams",
      {"SiTargetSourceIDs", "SiPadSourceIDs"},
      "Frame parameter name carrying source_id vector for each collection"};
  Gaudi::Property<std::vector<int>> m_detectorIDs{
      this, "DetectorIDs", {0, 1},
      "detector_id value (for logging) for each collection"};
  Gaudi::Property<std::vector<std::string>> m_ntupleNames{
      this, "NTupleNames", {"SiTarget", "SiPad"},
      "RNTuple name inside the output file for each collection"};
  Gaudi::Property<std::string> m_inputFile{
      this, "InputFile", "timewindows.root",
      "Input edm4hep file (same file read by IOSvc)"};
  Gaudi::Property<std::string> m_outputFile{
      this, "OutputFile", "ShipHits.root",
      "Output ROOT file containing all per-detector RNTuples"};

  // Path to the tracks.root file produced by ACTSProtoTracker.
  // If empty, track export is disabled.
  // The file must contain the same number of events as the hit input file.
  Gaudi::Property<std::string> m_trackFile{
      this, "TrackFile", "",
      "Path to tracks.root file from ACTSProtoTracker. "
      "Must have same number of events as hit input. "
      "Leave empty to skip track export."};

  Gaudi::Property<std::string> m_trackCollectionName{
    this, "TrackCollectionName", "ACTSTracks",
    "Name of the track collection inside TrackFile."};

  // Measurement collections (TrackerHit3D from converters)
  Gaudi::Property<std::vector<std::string>> m_measCollections{
      this, "MeasCollections", {},
      "TrackerHit3D collection names to export (e.g. SiTargetMeasurements)."};

  Gaudi::Property<std::vector<std::string>> m_measNtupleNames{
      this, "MeasNtupleNames", {},
      "RNTuple names for each MeasCollection (e.g. SiTargetMeas)."};

  Gaudi::Property<std::vector<std::string>> m_measBitFields{
      this, "MeasBitFields", {},
      "BitField strings for decoding cellID of each MeasCollection."};

  // One DataHandle per collection, built in initialize() from m_collections.
  mutable std::vector<
      std::unique_ptr<k4FWCore::DataHandle<edm4hep::SimCalorimeterHitCollection>>
  > m_inputHandles;

  // One BitFieldCoder per collection.
  mutable std::vector<
      std::unique_ptr<dd4hep::DDSegmentation::BitFieldCoder>
  > m_decoders;

  // Per-detector field names, populated in initialize() from the BitFieldCoders.
  mutable std::vector<std::vector<std::string>> m_detectorFieldNames;

  // podio::ROOTReader opened in initialize(), used in execute() to read
  // frame parameters (source_id vectors, t_window_start) by window index.
  mutable std::unique_ptr<podio::ROOTReader> m_reader;

  // Output TFile shared by all per-detector writers.
  mutable std::unique_ptr<TFile> m_outputTFile;

  // Per-detector RNTuple writers (one per collection).
  mutable std::vector<std::unique_ptr<ROOT::RNTupleWriter>> m_writers;

  // Per-detector field pointer bundles.
  struct DetectorFields {
    int*      fWindowID     = nullptr;
    uint64_t* fCellID       = nullptr;
    float*    fX            = nullptr;
    float*    fY            = nullptr;
    float*    fZ            = nullptr;
    float*    fE            = nullptr;
    float*    fT            = nullptr;
    int*      fSourceID     = nullptr;
    float*    fTWindowStart = nullptr;
    // Dynamic BitField field pointers for this detector only.
    std::map<std::string, int*>        fieldPtrs;
    // Keeps shared_ptrs alive after model is moved into writer.
    std::vector<std::shared_ptr<int>>  fieldSharedPtrs;
  };
  mutable std::vector<DetectorFields> m_detFields;

  // DataHandles for measurement collections (TrackerHit3D)
  mutable std::vector<
      std::unique_ptr<k4FWCore::DataHandle<edm4hep::TrackerHit3DCollection>>>
      m_measHandles;

  // RNTuple writers for measurements
  mutable std::vector<std::unique_ptr<ROOT::RNTupleWriter>> m_measWriters;

  struct MeasFields {
    int*      fWindowID = nullptr;
    float*    fX        = nullptr;
    float*    fY        = nullptr;
    float*    fZ        = nullptr;
    int*      fPlane    = nullptr;
    int*      fLayer    = nullptr;
    int*      fColumn   = nullptr;
    int*      fRow      = nullptr;
    float*    fEdep     = nullptr;
    int*      fSourceID = nullptr;
    std::vector<std::shared_ptr<int>>   intPtrs;
    std::vector<std::shared_ptr<float>> floatPtrs;
  };
  mutable std::vector<MeasFields> m_measFields;
  mutable std::vector<std::unique_ptr<dd4hep::DDSegmentation::BitFieldCoder>>
      m_measDecoders;

  // Direct podio reader for track file (bypasses Gaudi DataHandle)
  mutable std::unique_ptr<podio::ROOTReader> m_trackReader;
  mutable std::size_t                        m_nTrackEvents{0};

  // Track RNTuple field pointers (flat, no nested structs)
  mutable std::unique_ptr<ROOT::RNTupleWriter> m_trackWriter;
  mutable int*   m_fTrackWindowID = nullptr;
  mutable int*   m_fTrackID       = nullptr;
  mutable float* m_fTrackChi2     = nullptr;
  mutable int*   m_fTrackNdf      = nullptr;
  mutable float* m_fTrackSeedX    = nullptr;
  mutable float* m_fTrackSeedY    = nullptr;
  // Keep shared_ptrs alive after model is moved into writer
  mutable std::vector<std::shared_ptr<int>>   m_trackIntPtrs;
  mutable std::vector<std::shared_ptr<float>> m_trackFloatPtrs;

  // Per-surface track states RNTuple ("ACTSTrackStates"), one row per surface per track.
  mutable std::unique_ptr<ROOT::RNTupleWriter> m_trackStatesWriter;
  mutable int*   m_fTsWindowID  = nullptr;
  mutable int*   m_fTsTrackID   = nullptr;
  mutable int*   m_fTsStateIdx  = nullptr;
  mutable float* m_fTsX         = nullptr;
  mutable float* m_fTsY         = nullptr;
  mutable float* m_fTsZ         = nullptr;
  mutable float* m_fTsLoc0      = nullptr;
  mutable float* m_fTsTilt      = nullptr;
  mutable float* m_fTsChi2      = nullptr;
  mutable int*   m_fTsNmeas     = nullptr;
  mutable std::vector<std::shared_ptr<int>>   m_tsIntPtrs;
  mutable std::vector<std::shared_ptr<float>> m_tsFloatPtrs;

  // Window counter: incremented in each execute() call.
  mutable std::atomic<long long> m_windowCount{0};
};

DECLARE_COMPONENT(EDM4HEP2RNTuple)

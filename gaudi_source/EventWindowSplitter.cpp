#include "Gaudi/Algorithm.h"
#include "GaudiKernel/MsgStream.h"
#include "edm4hep/SimCalorimeterHitCollection.h"
#include "edm4hep/CaloHitContributionCollection.h"
#include "podio/ROOTReader.h"
#include "podio/ROOTWriter.h"
#include "podio/Frame.h"
#include "ValidationUtils.h"

#include <algorithm>
#include <limits>
#include <map>
#include <tuple>
#include <unordered_map>
#include <vector>

struct ContribData {
  float energy;
  float time;
  // source_id is encoded in the PDG field of CaloHitContribution by EventShuffler
  // because edm4hep::CaloHitContribution has no dedicated source field.
  int   source_id;
  int   pdg;         // real Geant4 PDG code, read from SiPadContribPDGs parameter
  // Exact Geant4 step position inside the sensitive volume
  float stepX, stepY, stepZ;
  // Centre of the parent SimCalorimeterHit sensitive volume
  float hitX, hitY, hitZ;
};

struct TimedContrib {
  uint64_t   cellID;
  ContribData data;
};

class EventWindowSplitter : public Gaudi::Algorithm {
public:
  EventWindowSplitter(const std::string& name, ISvcLocator* svcLoc)
      : Gaudi::Algorithm(name, svcLoc) {}

  StatusCode initialize() override {
    try {
      return Gaudi::Algorithm::initialize();
    } catch (const std::exception& e) {
      error() << "[EventWindowSplitter] Exception in initialize(): " << e.what() << endmsg;
      return StatusCode::FAILURE;
    } catch (...) {
      error() << "[EventWindowSplitter] Unknown exception in initialize()." << endmsg;
      return StatusCode::FAILURE;
    }
  }

  // execute() runs exactly once (EvtMax=1). All work is done here by reading
  // the input file directly via podio::ROOTReader (bypassing IOSvc).
  StatusCode execute(const EventContext&) const override {
    try {
      // --- 0. Open input super-event file ---
      podio::ROOTReader reader;
      reader.openFile(m_inputFile.value());

      if (reader.getEntries("events") == 0) {
        error() << "[EventWindowSplitter] Input file has 0 events." << endmsg;
        return StatusCode::FAILURE;
      }

      auto frame = podio::Frame(reader.readEntry("events", 0));

      // Validate collections and source ID parameters before processing
      bool ok = true;
      ok &= validateCollection(frame, m_inPadName.value(), "EventWindowSplitter", msgStream());
      if (!ok) return StatusCode::FAILURE;

      const auto& inPad = frame.get<edm4hep::SimCalorimeterHitCollection>(m_inPadName.value());

      validateParameter(frame, "SiPadSourceIDs", inPad.size(),
                        "EventWindowSplitter", msgStream());

      // Read source_id parameters written by EventShuffler
      std::vector<int> sourceIDsSiPad =
          frame.getParameter<std::vector<int>>("SiPadSourceIDs").value_or(std::vector<int>{});

      // Read per-contribution real Geant4 PDG codes written by EventShuffler.
      // value_or({}) provides backward compatibility with files that pre-date this feature.
      std::vector<int> pdgsSiPad =
          frame.getParameter<std::vector<int>>("SiPadContribPDGs").value_or(std::vector<int>{});

      auto makeSourceMap = [](const edm4hep::SimCalorimeterHitCollection& coll,
                              const std::vector<int>& ids) {
        std::unordered_map<uint64_t, int> m;
        size_t idx = 0;
        for (const auto& hit : coll)
          m[hit.getCellID()] = (idx < ids.size()) ? ids[idx++] : 0;
        return m;
      };

      // Build cellID -> source_id maps for fast lookup during window construction
      auto siPadSourceMap = makeSourceMap(inPad, sourceIDsSiPad);

      // --- 1. Build flat, sorted list of timed contributions ---
      auto buildList = [](const edm4hep::SimCalorimeterHitCollection& coll,
                          const std::vector<int>& pdgVec)
          -> std::vector<TimedContrib> {
        std::vector<TimedContrib> list;
        size_t pdgIdx = 0;
        for (const auto& hit : coll) {
          const auto& hp = hit.getPosition();
          for (const auto& contrib : hit.getContributions()) {
            const auto& sp = contrib.getStepPosition();
            const int realPDG = (pdgIdx < pdgVec.size()) ? pdgVec[pdgIdx] : 0;
            ++pdgIdx;
            list.push_back({hit.getCellID(),
                            {contrib.getEnergy(), contrib.getTime(),
                             contrib.getPDG(),  // source_id encoded in PDG field by EventShuffler
                             realPDG,
                             sp.x, sp.y, sp.z,
                             hp.x, hp.y, hp.z}});
          }
        }
        std::sort(list.begin(), list.end(),
                  [](const TimedContrib& a, const TimedContrib& b) {
                    return a.data.time < b.data.time;
                  });
        return list;
      };

      const auto siPadList = buildList(inPad, pdgsSiPad);

      if (siPadList.empty()) {
        error() << "[EventWindowSplitter] No contributions found in '"
                << m_inputFile.value() << "'. Cannot write time windows."
                << " Likely cause: simulation output has hits but no CaloHitContributions."
                << " Ensure the simulation was run WITHOUT SIM.part.userParticleHandler=''." << endmsg;
        return StatusCode::FAILURE;
      }

      // --- 2. Determine global t_start ---
      float t_start     = siPadList.front().data.time;
      float t_end_global = siPadList.back().data.time;

      // --- 3. Assign each contribution to a window index ---
      // Window k covers [t_start + k*W, t_start + (k+1)*W).
      const float W = static_cast<float>(m_windowSize);

      auto distribute =
          [&](const std::vector<TimedContrib>& list)
          -> std::map<int, std::map<uint64_t, std::vector<ContribData>>> {
        std::map<int, std::map<uint64_t, std::vector<ContribData>>> windows;
        for (const auto& tc : list) {
          const int wIdx = std::max(0, static_cast<int>((tc.data.time - t_start) / W));
          windows[wIdx][tc.cellID].push_back(tc.data);
        }
        return windows;
      };

      const auto windowsSiPad = distribute(siPadList);

      int numWindows = 0;
      if (!windowsSiPad.empty()) numWindows = windowsSiPad.rbegin()->first + 1;

      // --- 4. Write one frame per window ---
      podio::ROOTWriter writer(m_outputFile.value());
      long long totalContribs = 0;

      // Returns hits collection, contributions collection, parallel source_id vector,
      // and parallel per-contribution real PDG vector.
      auto buildWindowColl =
          [&](int wIdx,
              float t_window_start,
              const std::map<int, std::map<uint64_t, std::vector<ContribData>>>& windows,
              const std::unordered_map<uint64_t, int>& sourceMap)
          -> std::tuple<edm4hep::SimCalorimeterHitCollection,
                        edm4hep::CaloHitContributionCollection,
                        std::vector<int>,
                        std::vector<int>> {
        edm4hep::SimCalorimeterHitCollection   hits;
        edm4hep::CaloHitContributionCollection contribs;
        std::vector<int> sourceIDs;
        std::vector<int> contribPDGs;

        auto it = windows.find(wIdx);
        if (it == windows.end()) return {std::move(hits), std::move(contribs), std::move(sourceIDs), std::move(contribPDGs)};

        for (const auto& [cellID, contribVec] : it->second) {
          auto hit = hits.create();
          hit.setCellID(cellID);
          hit.setPosition({contribVec[0].hitX, contribVec[0].hitY, contribVec[0].hitZ});

          const int hitSourceID = sourceMap.count(cellID) ? sourceMap.at(cellID) : 0;
          sourceIDs.push_back(hitSourceID);

          float totalE = 0.0f;
          for (const auto& cd : contribVec) {
            auto contrib = contribs.create();
            contrib.setEnergy(cd.energy);
            contrib.setTime(cd.time - t_window_start);
            contrib.setPDG(cd.source_id);
            contrib.setStepPosition({cd.stepX, cd.stepY, cd.stepZ});
            hit.addToContributions(contrib);
            contribPDGs.push_back(cd.pdg);
            totalE += cd.energy;
            ++totalContribs;
          }
          hit.setEnergy(totalE);
        }
        return {std::move(hits), std::move(contribs), std::move(sourceIDs), std::move(contribPDGs)};
      };

      for (int wIdx = 0; wIdx < numWindows; ++wIdx) {
        const float t_window_start =
            t_start + static_cast<float>(wIdx) * static_cast<float>(m_windowSize);
        podio::Frame windowFrame;

        auto [siPadColl, siPadContribs, wSrcPad, wPdgPad] =
            buildWindowColl(wIdx, t_window_start, windowsSiPad, siPadSourceMap);

        windowFrame.put(std::move(siPadColl),    m_outPadName.value());
        windowFrame.put(std::move(siPadContribs), m_outPadName.value() + "_Contributions");

        windowFrame.putParameter("SiPadSourceIDs",    wSrcPad);
        windowFrame.putParameter("SiPadContribPDGs",  wPdgPad);
        windowFrame.putParameter("t_window_start",    t_window_start);

        writer.writeFrame(windowFrame, "events");
      }
      writer.finish();

      // --- 5. Summary ---
      const double avgPerWindow = (numWindows > 0)
          ? static_cast<double>(siPadList.size()) / numWindows : 0.0;

      info() << "[EventWindowSplitter] Total windows written: " << numWindows << endmsg;
      info() << "[EventWindowSplitter] Time range: " << t_start << " to " << t_end_global << " ns" << endmsg;
      info() << "[EventWindowSplitter] Avg contributions per window: " << avgPerWindow << endmsg;

      return StatusCode::SUCCESS;
    } catch (const std::exception& e) {
      error() << "[EventWindowSplitter] Exception in execute(): " << e.what() << endmsg;
      return StatusCode::FAILURE;
    } catch (...) {
      error() << "[EventWindowSplitter] Unknown exception in execute()." << endmsg;
      return StatusCode::FAILURE;
    }
  }

  StatusCode finalize() override {
    try {
      return Gaudi::Algorithm::finalize();
    } catch (const std::exception& e) {
      error() << "[EventWindowSplitter] Exception in finalize(): " << e.what() << endmsg;
      return StatusCode::FAILURE;
    } catch (...) {
      error() << "[EventWindowSplitter] Unknown exception in finalize()." << endmsg;
      return StatusCode::FAILURE;
    }
  }

private:
  Gaudi::Property<std::string> m_inputFile{
      this, "InputFile", "shuffled.root", "Input super-event file"};
  Gaudi::Property<std::string> m_inPadName{
      this, "InputCollectionSiPad", "SiPadHitsMerged", "Input SiPad collection"};
  Gaudi::Property<std::string> m_outputFile{
      this, "OutputFile", "timewindows.edm4hep.root", "Output ROOT file for time-windowed events"};
  Gaudi::Property<std::string> m_outPadName{
      this, "OutputCollectionSiPad", "SiPadHitsWindowed", "Output SiPad collection name"};
  Gaudi::Property<double> m_windowSize{
      this, "WindowSize", 25.0, "Time window size (ns)"};
};

DECLARE_COMPONENT(EventWindowSplitter)

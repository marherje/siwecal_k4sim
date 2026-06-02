#include "Gaudi/Algorithm.h"
#include "GaudiKernel/MsgStream.h"
#include "k4FWCore/DataHandle.h"
#include "edm4hep/SimCalorimeterHitCollection.h"
#include "edm4hep/CaloHitContributionCollection.h"
#include "podio/ROOTReader.h"
#include "podio/ROOTWriter.h"
#include "podio/Frame.h"
#include "ValidationUtils.h"

#include <algorithm>
#include <map>
#include <numeric>
#include <string>
#include <tuple>
#include <vector>

class EventShuffler : public Gaudi::Algorithm {
public:
  EventShuffler(const std::string& name, ISvcLocator* svcLoc)
      : Gaudi::Algorithm(name, svcLoc) {}

  StatusCode initialize() override {
    try {
      StatusCode sc = Gaudi::Algorithm::initialize();
      if (sc.isFailure()) return sc;

      const auto& files  = m_inputFiles.value();
      const auto& ids    = m_sourceIDs.value();
      const auto& delays = m_delays.value();
      const auto& colsP  = m_collsSiPad.value();

      if (files.size() != ids.size()   ||
          files.size() != delays.size() ||
          files.size() != colsP.size()) {
        error() << "[EventShuffler] InputFiles, SourceIDs, Delays, CollectionsSiPad"
                << " must all have the same size. Got: "
                << files.size() << ", " << ids.size() << ", "
                << delays.size() << ", " << colsP.size()
                << endmsg;
        return StatusCode::FAILURE;
      }

      info() << "[EventShuffler] Configured with " << files.size() << " input sources:" << endmsg;
      for (size_t i = 0; i < files.size(); ++i) {
        info() << "[EventShuffler]   source_id=" << ids[i]
               << "  delay=" << delays[i] << " ns"
               << "  file=" << files[i]
               << endmsg;
      }
      return sc;
    } catch (const std::exception& e) {
      error() << "[EventShuffler] Exception in initialize(): " << e.what() << endmsg;
      return StatusCode::FAILURE;
    } catch (...) {
      error() << "[EventShuffler] Unknown exception in initialize()." << endmsg;
      return StatusCode::FAILURE;
    }
  }

  // execute() is a no-op: all work is done in finalize() by reading the
  // input files directly via podio::ROOTReader (bypassing IOSvc).
  StatusCode execute(const EventContext&) const override {
    return StatusCode::SUCCESS;
  }

  StatusCode finalize() override {
    try {
      const auto& files  = m_inputFiles.value();
      const auto& ids    = m_sourceIDs.value();
      const auto& delays = m_delays.value();
      const auto& colsP  = m_collsSiPad.value();

      // Quick scan to get event counts per source (needed for exhaustion warning).
      std::vector<size_t> nEventsVec(files.size());
      for (size_t i = 0; i < files.size(); ++i) {
        podio::ROOTReader r;
        r.openFile(files[i]);
        nEventsVec[i] = r.getEntries("events");
      }
      const size_t maxNEvents = *std::max_element(nEventsVec.begin(), nEventsVec.end());

      // Main processing: validate each source, then accumulate contributions.
      for (size_t i = 0; i < files.size(); ++i) {
        const size_t nEvents  = nEventsVec[i];
        const int    sourceID = ids[i];
        const double delay    = delays[i];

        if (nEvents == 0) {
          warning() << "[EventShuffler] Source source_id=" << sourceID
                    << " file=" << files[i] << " has 0 events. Skipping." << endmsg;
          continue;
        }

        info() << "[EventShuffler] Source source_id=" << sourceID
               << ": file=" << files[i] << "  events=" << nEvents << endmsg;

        podio::ROOTReader reader;
        reader.openFile(files[i]);

        // Validate first frame before processing all events
        {
          auto firstFrame = podio::Frame(reader.readEntry("events", 0));
          bool ok = true;
          ok &= validateCollection(firstFrame, colsP[i], "EventShuffler", msgStream());
          if (!ok) {
            error() << "[EventShuffler] Validation failed for file: "
                    << files[i] << ". Aborting." << endmsg;
            return StatusCode::FAILURE;
          }
          info() << "[EventShuffler] Validation passed for file: "
                 << files[i] << "  events=" << nEvents << endmsg;
        }

        // Process all events (readEntry uses absolute index, so re-reading from 0 is fine)
        const size_t maxEv = (m_maxEventsPerSource.value() > 0)
                             ? std::min(nEvents, static_cast<size_t>(m_maxEventsPerSource.value()))
                             : nEvents;
        if (maxEv < nEvents) {
          info() << "[EventShuffler] source_id=" << sourceID
                 << ": limiting to " << maxEv << " of " << nEvents << " events." << endmsg;
        }
        for (size_t ev = 0; ev < maxEv; ++ev) {
          auto frame = podio::Frame(reader.readEntry("events", ev));

          const auto& hitsPad =
              frame.get<edm4hep::SimCalorimeterHitCollection>(colsP[i]);

          auto accum = [&](const edm4hep::SimCalorimeterHitCollection& hits,
                           std::map<uint64_t, std::vector<ContribData>>& buf) {
            for (const auto& hit : hits) {
              auto& vec = buf[hit.getCellID()];
              const auto& hp = hit.getPosition();
              for (const auto& contrib : hit.getContributions()) {
                const auto& sp = contrib.getStepPosition();
                ContribData cd;
                cd.energy    = contrib.getEnergy();
                cd.time      = static_cast<float>(static_cast<double>(ev) * delay)
                               + contrib.getTime();
                // ddsim leaves CaloHitContribution.PDG = 0; pull it from the linked MCParticle instead.
                const auto mcp = contrib.getParticle();
                cd.pdg       = mcp.isAvailable() ? mcp.getPDG() : 0;
                cd.source_id = sourceID;
                cd.stepX     = sp.x;
                cd.stepY     = sp.y;
                cd.stepZ     = sp.z;
                cd.hitX      = hp.x;
                cd.hitY      = hp.y;
                cd.hitZ      = hp.z;
                vec.push_back(cd);
              }
            }
          };

          accum(hitsPad, m_bufferSiPad);
        }

        if (nEvents < maxNEvents) {
          warning() << "[EventShuffler] WARNING: source_id=" << sourceID
                    << " exhausted after " << nEvents << " events. "
                    << "Other sources may have more events. "
                    << "Check your simulation inputs." << endmsg;
        }
      }

      // Build the output super-event.
      auto buildColl =
          [](const std::map<uint64_t, std::vector<ContribData>>& buf)
          -> std::tuple<edm4hep::SimCalorimeterHitCollection,
                        edm4hep::CaloHitContributionCollection,
                        std::vector<int>,
                        std::vector<int>> {
        edm4hep::SimCalorimeterHitCollection   hits;
        edm4hep::CaloHitContributionCollection contribs;
        std::vector<int> sourceIDs;
        std::vector<int> pdgs;  // per-contribution real Geant4 PDG codes
        for (const auto& [cellID, vec] : buf) {
          auto outHit = hits.create();
          outHit.setCellID(cellID);
          outHit.setPosition({vec[0].hitX, vec[0].hitY, vec[0].hitZ});

          // Determine hit-level source_id: common value if uniform, else 0 (mixed).
          int hitSourceID = vec[0].source_id;
          for (const auto& cd : vec) {
            if (cd.source_id != hitSourceID) { hitSourceID = 0; break; }
          }
          sourceIDs.push_back(hitSourceID);

          float totalE = 0.0f;
          for (const auto& cd : vec) {
            // source_id is encoded in the PDG field because
            // edm4hep::CaloHitContribution has no dedicated source field.
            // PDG is not used downstream in this pipeline.
            auto contrib = contribs.create();
            contrib.setEnergy(cd.energy);
            contrib.setTime(cd.time);
            contrib.setPDG(cd.source_id);
            contrib.setStepPosition({cd.stepX, cd.stepY, cd.stepZ});
            outHit.addToContributions(contrib);
            pdgs.push_back(cd.pdg);  // parallel PDG vector (per contribution)
            totalE += cd.energy;
          }
          outHit.setEnergy(totalE);
        }
        return {std::move(hits), std::move(contribs), std::move(sourceIDs), std::move(pdgs)};
      };

      auto [outSiPad, outSiPadContribs, sourceIDsSiPad, pdgsSiPad] = buildColl(m_bufferSiPad);

      podio::ROOTWriter writer(m_outputFile.value());
      {
        podio::Frame outFrame;
        outFrame.put(std::move(outSiPad),         m_outputCollSiPad.value());
        outFrame.put(std::move(outSiPadContribs), m_outputCollSiPad.value() + "_Contributions");
        outFrame.putParameter("SiPadSourceIDs",    sourceIDsSiPad);
        // Per-contribution real Geant4 PDG codes (parallel to *_Contributions collections)
        outFrame.putParameter("SiPadContribPDGs",  pdgsSiPad);
        writer.writeFrame(outFrame, "events");
      }
      writer.finish();

      info() << "[EventShuffler] Super-event written to: " << m_outputFile.value() << endmsg;
      info() << "[EventShuffler] SiPad cellIDs: " << outSiPad.size() << endmsg;

      return Gaudi::Algorithm::finalize();
    } catch (const std::exception& e) {
      error() << "[EventShuffler] Exception in finalize(): " << e.what() << endmsg;
      return StatusCode::FAILURE;
    } catch (...) {
      error() << "[EventShuffler] Unknown exception in finalize()." << endmsg;
      return StatusCode::FAILURE;
    }
  }

private:
  struct ContribData {
    float energy;
    float time;        // absolute time in ns
    int   source_id;
    int   pdg;         // real Geant4 PDG code from the input ddsim file
    float stepX, stepY, stepZ;
    float hitX, hitY, hitZ;
  };

  // Python-configurable properties
  Gaudi::Property<std::vector<std::string>> m_inputFiles{
      this, "InputFiles", {}, "List of input edm4hep ROOT files"};
  Gaudi::Property<std::vector<int>> m_sourceIDs{
      this, "SourceIDs", {}, "Source ID for each input file (same order as InputFiles)"};
  Gaudi::Property<std::vector<double>> m_delays{
      this, "Delays", {}, "Inter-event delay in ns for each input file"};
  Gaudi::Property<std::vector<std::string>> m_collsSiPad{
      this, "CollectionsSiPad", {}, "SiPad collection name for each input file"};
  Gaudi::Property<std::string> m_outputFile{
      this, "OutputFile", "shuffled.root", "Output super-event ROOT file"};
  Gaudi::Property<std::string> m_outputCollSiPad{
      this, "OutputCollectionSiPad", "SiPadHitsMerged", "Output SiPad collection name"};
  Gaudi::Property<int> m_maxEventsPerSource{
      this, "MaxEventsPerSource", 0,
      "Maximum number of events to read from each source (0 = no limit)"};

  mutable std::map<uint64_t, std::vector<ContribData>> m_bufferSiPad;
};

DECLARE_COMPONENT(EventShuffler)

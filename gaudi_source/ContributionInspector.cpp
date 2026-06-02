#include "ContributionInspector.h"

#include <limits>

ContributionInspector::ContributionInspector(const std::string& name, ISvcLocator* svcLoc)
    : Gaudi::Algorithm(name, svcLoc) {}

StatusCode ContributionInspector::initialize() {
  try {
    StatusCode sc = Gaudi::Algorithm::initialize();
    if (sc.isFailure()) return sc;
    m_targetHandle = std::make_unique<k4FWCore::DataHandle<edm4hep::SimCalorimeterHitCollection>>(
        m_targetName.value(), Gaudi::DataHandle::Reader, this);
    m_pixelHandle  = std::make_unique<k4FWCore::DataHandle<edm4hep::SimCalorimeterHitCollection>>(
        m_pixelName.value(),  Gaudi::DataHandle::Reader, this);
    return sc;
  } catch (const std::exception& e) {
    error() << "[ContributionInspector] Exception in initialize(): " << e.what() << endmsg;
    return StatusCode::FAILURE;
  } catch (...) {
    error() << "[ContributionInspector] Unknown exception in initialize()." << endmsg;
    return StatusCode::FAILURE;
  }
}

StatusCode ContributionInspector::execute(const EventContext&) const {
  try {
    inspectCollection(*m_targetHandle->get(), m_targetName.value());
    inspectCollection(*m_pixelHandle->get(),  m_pixelName.value());
    return StatusCode::SUCCESS;
  } catch (const std::exception& e) {
    error() << "[ContributionInspector] Exception in execute(): " << e.what() << endmsg;
    return StatusCode::FAILURE;
  } catch (...) {
    error() << "[ContributionInspector] Unknown exception in execute()." << endmsg;
    return StatusCode::FAILURE;
  }
}

StatusCode ContributionInspector::finalize() {
  try {
    m_targetHandle.reset();
    m_pixelHandle.reset();
    return Gaudi::Algorithm::finalize();
  } catch (const std::exception& e) {
    error() << "[ContributionInspector] Exception in finalize(): " << e.what() << endmsg;
    return StatusCode::FAILURE;
  } catch (...) {
    error() << "[ContributionInspector] Unknown exception in finalize()." << endmsg;
    return StatusCode::FAILURE;
  }
}

void ContributionInspector::inspectCollection(
    const edm4hep::SimCalorimeterHitCollection& hits,
    const std::string& collName) const
{
  info() << "[" << collName << "] Event has " << hits.size() << " hits" << endmsg;

  int   totalContribs    = 0;
  int   zeroTimeContribs = 0;
  int   posTimeContribs  = 0;
  float tMin = std::numeric_limits<float>::max();
  float tMax = std::numeric_limits<float>::lowest();

  for (int i = 0; i < static_cast<int>(hits.size()); ++i) {
    const auto& hit       = hits[i];
    const bool  doPrint   = (i < m_maxHitsToPrint);
    int         j         = 0;

    for (const auto& contrib : hit.getContributions()) {
      const float t   = contrib.getTime();
      const float e   = contrib.getEnergy();
      const int   pdg = contrib.getPDG();

      ++totalContribs;
      if (t == 0.0f) ++zeroTimeContribs;
      else if (t > 0.0f) ++posTimeContribs;
      if (t < tMin) tMin = t;
      if (t > tMax) tMax = t;

      if (doPrint) {
        debug() << "[" << collName << "]"
                << " hit=" << i << " contrib=" << j
                << " time=" << t << " ns"
                << " energy=" << e << " GeV"
                << " pdg=" << pdg
                << endmsg;
      }
      ++j;
    }
  }

  // Per-collection summary
  info() << "[" << collName << "] Total contributions: " << totalContribs << endmsg;

  if (totalContribs > 0) {
    const double zeroFrac = 100.0 * zeroTimeContribs / totalContribs;
    const double posFrac  = 100.0 * posTimeContribs  / totalContribs;
    info() << "[" << collName << "] Time range: " << tMin << " ns to " << tMax << " ns" << endmsg;
    info() << "[" << collName << "] Contributions with time=0: "
           << zeroTimeContribs << " (" << zeroFrac << "%)" << endmsg;
    info() << "[" << collName << "] Contributions with time>0: "
           << posTimeContribs  << " (" << posFrac  << "%)" << endmsg;

    if (zeroTimeContribs > totalContribs / 2) {
      warning() << "[" << collName << "] More than 50% of contributions have time=0 -- "
                << "the time field may not be populated; "
                << "time-based event merger cannot be implemented." << endmsg;
    }
  } else {
    info() << "[" << collName << "] Time range: N/A (no contributions)" << endmsg;
    info() << "[" << collName << "] Contributions with time=0: 0 (N/A)" << endmsg;
    info() << "[" << collName << "] Contributions with time>0: 0 (N/A)" << endmsg;
  }
}

DECLARE_COMPONENT(ContributionInspector)

#pragma once

#include "Gaudi/Algorithm.h"
#include "GaudiKernel/MsgStream.h"
#include "k4FWCore/DataHandle.h"
#include "edm4hep/SimCalorimeterHitCollection.h"

#include <memory>

// Diagnostic algorithm that inspects the time structure of CaloHitContributions.
// No output collection is written. Intended to verify that the Geant4 time-of-flight
// field is populated before implementing a time-based event merger.
class ContributionInspector : public Gaudi::Algorithm {
public:
  ContributionInspector(const std::string& name, ISvcLocator* svcLoc);

  StatusCode initialize() override;
  StatusCode execute(const EventContext&) const override;
  StatusCode finalize() override;

private:
  void inspectCollection(const edm4hep::SimCalorimeterHitCollection& hits,
                         const std::string& collName) const;

  Gaudi::Property<std::string> m_pixelName{
      this, "SiPadCollection", "SiPadHits", "Input SiPad SimCalorimeterHit collection"};
  Gaudi::Property<int> m_maxHitsToPrint{
      this, "MaxHitsToPrint", 20, "Max number of hits per collection for which contributions are printed at DEBUG level"};

  mutable std::unique_ptr<k4FWCore::DataHandle<edm4hep::SimCalorimeterHitCollection>> m_pixelHandle;
};

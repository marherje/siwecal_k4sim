#include "DDG4/Factories.h"
#include "DDG4/Geant4SteppingAction.h"

#include "G4LogicalVolume.hh"
#include "G4ParticleDefinition.hh"
#include "G4Step.hh"
#include "G4SystemOfUnits.hh"
#include "G4Track.hh"
#include "G4UserLimits.hh"
#include "G4VPhysicalVolume.hh"

#include <algorithm>
#include <cstdint>
#include <iostream>
#include <limits>
#include <string>

namespace dd4hep::sim {

class StepLengthDebugger : public Geant4SteppingAction {
public:
  StepLengthDebugger(Geant4Context* context, const std::string& name)
      : Geant4SteppingAction(context, name) {
    declareProperty("VolumeSubstring", m_volumeSubstring);
    declareProperty("MaxPrintedSteps", m_maxPrintedSteps);
    declareProperty("ToleranceMm", m_toleranceMm);
    declareProperty("PrintOnlyOverLimit", m_printOnlyOverLimit);
  }

  ~StepLengthDebugger() override {
    std::cout << "[StepLengthDebugger] summary: total_steps=" << m_totalSteps
              << " matched_steps=" << m_matchedSteps
              << " with_limits=" << m_withLimits
              << " without_limits=" << m_withoutLimits
              << " over_limit=" << m_overLimit
              << " max_step=" << (m_maxStep / mm) << " mm"
              << " max_limited_step=" << (m_maxLimitedStep / mm) << " mm"
              << " max_allowed=" << (m_maxAllowed / mm) << " mm"
              << std::endl;
  }

  void operator()(const G4Step* step, G4SteppingManager*) override {
    ++m_totalSteps;
    if (step == nullptr) {
      return;
    }

    const auto* preStepPoint = step->GetPreStepPoint();
    if (preStepPoint == nullptr) {
      return;
    }

    auto* physicalVolume = preStepPoint->GetPhysicalVolume();
    if (physicalVolume == nullptr) {
      return;
    }

    auto* logicalVolume = physicalVolume->GetLogicalVolume();
    if (logicalVolume == nullptr) {
      return;
    }

    const std::string physicalName = physicalVolume->GetName();
    const std::string logicalName = logicalVolume->GetName();
    if (!matchesVolume(physicalName) && !matchesVolume(logicalName)) {
      return;
    }

    ++m_matchedSteps;

    const G4double stepLength = step->GetStepLength();
    m_maxStep = std::max(m_maxStep, stepLength);

    auto* limits = logicalVolume->GetUserLimits();
    G4double maxAllowedStep = std::numeric_limits<G4double>::infinity();
    bool hasLimits = false;
    bool overLimit = false;

    if (limits != nullptr) {
      hasLimits = true;
      ++m_withLimits;
      m_maxLimitedStep = std::max(m_maxLimitedStep, stepLength);
      maxAllowedStep = limits->GetMaxAllowedStep(*step->GetTrack());
      m_maxAllowed = std::max(m_maxAllowed, maxAllowedStep);
      overLimit = stepLength > maxAllowedStep + m_toleranceMm * mm;
      if (overLimit) {
        ++m_overLimit;
      }
    } else {
      ++m_withoutLimits;
    }

    const bool shouldPrint =
        m_printedSteps < m_maxPrintedSteps &&
        (!m_printOnlyOverLimit || overLimit || !hasLimits);
    if (!shouldPrint) {
      return;
    }

    ++m_printedSteps;
    const auto* track = step->GetTrack();
    const auto* particle = track ? track->GetParticleDefinition() : nullptr;
    const std::string particleName = particle ? particle->GetParticleName() : "";
    std::cout << "[StepLengthDebugger] step=" << m_matchedSteps
              << " particle=" << particleName
              << " physical=" << physicalName
              << " logical=" << logicalName
              << " step_length=" << (stepLength / mm) << " mm"
              << " max_allowed=" << (maxAllowedStep / mm) << " mm"
              << " has_limits=" << (hasLimits ? 1 : 0)
              << " over_limit=" << (overLimit ? 1 : 0)
              << std::endl;
  }

private:
  bool matchesVolume(const std::string& name) const {
    return m_volumeSubstring.empty() ||
           name.find(m_volumeSubstring) != std::string::npos;
  }

  std::string m_volumeSubstring{"SiPad"};
  int m_maxPrintedSteps{40};
  double m_toleranceMm{1e-6};
  bool m_printOnlyOverLimit{false};

  std::int64_t m_totalSteps{0};
  std::int64_t m_matchedSteps{0};
  std::int64_t m_printedSteps{0};
  std::int64_t m_withLimits{0};
  std::int64_t m_withoutLimits{0};
  std::int64_t m_overLimit{0};
  G4double m_maxStep{0.0};
  G4double m_maxLimitedStep{0.0};
  G4double m_maxAllowed{0.0};
};

} // namespace dd4hep::sim

DECLARE_GEANT4ACTION(StepLengthDebugger)

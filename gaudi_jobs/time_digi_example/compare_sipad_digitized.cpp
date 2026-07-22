#include "CellShaping.hh"

#include "edm4hep/SimCalorimeterHitCollection.h"
#include "podio/Frame.h"
#include "podio/ROOTReader.h"
#include "podio/UserDataCollection.h"

#include <cmath>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <limits>
#include <numeric>
#include <random>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

constexpr std::uint64_t kEventSeedStride = 0x9e3779b97f4a7c15ULL;

void printUsage(const char* argv0) {
  std::cout
      << "Usage: " << argv0
      << " [input.root] [digitized.root] [collection] [event] [max_hits] [options]\n"
      << "\nOptions:\n"
      << "  --digi-collection NAME    Digitized hit collection (default: SiPadHitsDigi)\n"
      << "  --digi-energy NAME        Digitized energy UserData collection\n"
      << "                            (default: SiPadHitsDigiDigitizedEnergy)\n"
      << "  --digi-time NAME          Digitized time UserData collection\n"
      << "                            (default: SiPadHitsDigiDigitizedTime)\n"
      << "  --mip-value GEV           MIP calibration [GeV/MIP] (default: 0.0002)\n"
      << "  --threshold MIP           Shaping threshold [MIP] (default: 0.5)\n"
      << "  --delay NS                Slow sample delay after fast trigger [ns] (default: 160)\n"
      << "  --tau-fast NS             Fast shaping time [ns] (default: 30)\n"
      << "  --tau-slow NS             Slow shaping time [ns] (default: 180)\n"
      << "  --order-fast N            Fast CR-RC order (default: 2)\n"
      << "  --order-slow N            Slow CR-RC order (default: 2)\n"
      << "  --fast-window NS          Fast peak search window [ns] (default: 200)\n"
      << "  --slow-window NS          Slow peak search window [ns] (default: 500)\n"
      << "  --fast-noise MIP          Fast noise sigma [MIP] (default: 1/30)\n"
      << "  --slow-noise MIP          Slow noise sigma [MIP] (default: 1/12)\n"
      << "  --peak-bins N             Coarse bins for peak search (default: 64)\n"
      << "  --refine-iterations N     Peak/threshold refinement iterations (default: 48)\n"
      << "  --trigger-bins N          Coarse bins for trigger search (default: 64)\n"
      << "  --seed N                  Base random seed (default: 5489)\n"
      << "  --help                    Show this message\n";
}

double parseDouble(const std::string& value, const std::string& option) {
  try {
    return std::stod(value);
  } catch (const std::exception&) {
    throw std::invalid_argument("Invalid value for " + option + ": " + value);
  }
}

int parseInt(const std::string& value, const std::string& option) {
  try {
    return std::stoi(value);
  } catch (const std::exception&) {
    throw std::invalid_argument("Invalid value for " + option + ": " + value);
  }
}

std::uint64_t parseUInt64(const std::string& value, const std::string& option) {
  try {
    return std::stoull(value);
  } catch (const std::exception&) {
    throw std::invalid_argument("Invalid value for " + option + ": " + value);
  }
}

std::string requireValue(int& i, int argc, char** argv, const std::string& option) {
  if (i + 1 >= argc) {
    throw std::invalid_argument("Missing value after " + option);
  }
  return argv[++i];
}

struct CompareConfig {
  std::string digiCollection = "SiPadHitsDigi";
  std::string digiEnergyCollection = "SiPadHitsDigiDigitizedEnergy";
  std::string digiTimeCollection = "SiPadHitsDigiDigitizedTime";
  std::uint64_t randomSeed = 5489ULL;
  siwecal::CellShapingConfig shaping;
};

struct HitSteps {
  std::size_t hitIndex = 0;
  std::uint64_t cellID = 0;
  double hitEnergyGeV = 0.0;
  std::vector<double> energyGeV;
  std::vector<double> timeNs;
};

HitSteps collectHitSteps(const edm4hep::SimCalorimeterHitCollection& hits,
                         std::size_t hitIndex) {
  const auto& hit = hits[hitIndex];
  HitSteps steps;
  steps.hitIndex = hitIndex;
  steps.cellID = hit.getCellID();
  steps.hitEnergyGeV = hit.getEnergy();
  steps.energyGeV.reserve(hit.contributions_size());
  steps.timeNs.reserve(hit.contributions_size());

  for (const auto& contrib : hit.getContributions()) {
    steps.energyGeV.push_back(contrib.getEnergy());
    steps.timeNs.push_back(contrib.getTime());
  }
  if (steps.energyGeV.empty()) {
    steps.energyGeV.push_back(hit.getEnergy());
    steps.timeNs.push_back(0.0);
  }
  return steps;
}

siwecal::FastSearchCellShapingResult runFastSearch(const HitSteps& steps,
                                                   const CompareConfig& cfg,
                                                   std::size_t eventIndex) {
  const auto seed = cfg.randomSeed +
                    kEventSeedStride * static_cast<std::uint64_t>(eventIndex + 1) +
                    static_cast<std::uint64_t>(steps.hitIndex);
  std::mt19937_64 rng(seed);
  return siwecal::fastSearchCellSteps(steps.energyGeV, steps.timeNs,
                                      cfg.shaping, rng);
}

bool passes(const siwecal::FastSearchCellShapingResult& result,
            const siwecal::CellShapingConfig& cfg) {
  return result.triggerTime >= 0.0 &&
         result.slowSignalSample >= cfg.mipThreshold;
}

void parseOptions(int argc, char** argv, CompareConfig& cfg) {
  cfg.shaping.mipValueGeV = 0.0002;
  for (int i = 6; i < argc; ++i) {
    const std::string opt = argv[i];
    if (opt == "--digi-collection") {
      cfg.digiCollection = requireValue(i, argc, argv, opt);
    } else if (opt == "--digi-energy") {
      cfg.digiEnergyCollection = requireValue(i, argc, argv, opt);
    } else if (opt == "--digi-time") {
      cfg.digiTimeCollection = requireValue(i, argc, argv, opt);
    } else if (opt == "--mip-value") {
      cfg.shaping.mipValueGeV = parseDouble(requireValue(i, argc, argv, opt), opt);
    } else if (opt == "--threshold") {
      cfg.shaping.mipThreshold = parseDouble(requireValue(i, argc, argv, opt), opt);
    } else if (opt == "--delay") {
      cfg.shaping.delayNs = parseDouble(requireValue(i, argc, argv, opt), opt);
    } else if (opt == "--tau-fast") {
      cfg.shaping.tauFastNs = parseDouble(requireValue(i, argc, argv, opt), opt);
    } else if (opt == "--tau-slow") {
      cfg.shaping.tauSlowNs = parseDouble(requireValue(i, argc, argv, opt), opt);
    } else if (opt == "--order-fast") {
      cfg.shaping.orderFast = parseInt(requireValue(i, argc, argv, opt), opt);
    } else if (opt == "--order-slow") {
      cfg.shaping.orderSlow = parseInt(requireValue(i, argc, argv, opt), opt);
    } else if (opt == "--fast-window") {
      cfg.shaping.fastWindowNs = parseDouble(requireValue(i, argc, argv, opt), opt);
    } else if (opt == "--slow-window") {
      cfg.shaping.slowWindowNs = parseDouble(requireValue(i, argc, argv, opt), opt);
    } else if (opt == "--fast-noise") {
      cfg.shaping.fastNoiseMIP = parseDouble(requireValue(i, argc, argv, opt), opt);
    } else if (opt == "--slow-noise") {
      cfg.shaping.slowNoiseMIP = parseDouble(requireValue(i, argc, argv, opt), opt);
    } else if (opt == "--peak-bins") {
      cfg.shaping.peakSearchBins = parseInt(requireValue(i, argc, argv, opt), opt);
    } else if (opt == "--refine-iterations") {
      cfg.shaping.refineIterations = parseInt(requireValue(i, argc, argv, opt), opt);
    } else if (opt == "--trigger-bins") {
      cfg.shaping.triggerSearchBins = parseInt(requireValue(i, argc, argv, opt), opt);
    } else if (opt == "--seed") {
      cfg.randomSeed = parseUInt64(requireValue(i, argc, argv, opt), opt);
    } else {
      throw std::invalid_argument("Unknown option: " + opt);
    }
  }
}

} // namespace

int main(int argc, char** argv) {
  for (int i = 1; i < argc; ++i) {
    if (std::string(argv[i]) == "--help") {
      printUsage(argv[0]);
      return 0;
    }
  }

  const std::string inputFile =
      argc > 1 ? argv[1] :
      "/home/llr/ilc/shi/data/siwecal_k4sim/output/output_PG_gamma_10GeV_100evt.edm4hep.root";
  const std::string digitizedFile =
      argc > 2 ? argv[2] :
      "/home/llr/ilc/shi/data/siwecal_k4sim/output/output_PG_gamma_10GeV_100evt_real_digitized.edm4hep.root";
  const std::string inputCollection = argc > 3 ? argv[3] : "SiPadHits";
  const std::size_t eventIndex = argc > 4 ? std::stoul(argv[4]) : 0;
  const int maxHits = argc > 5 ? std::stoi(argv[5]) : 10;

  CompareConfig cfg;
  parseOptions(argc, argv, cfg);

  podio::ROOTReader inputReader;
  inputReader.openFile(inputFile);
  if (eventIndex >= inputReader.getEntries("events")) {
    throw std::runtime_error("Requested event index is outside input file");
  }
  auto inputFrame = podio::Frame(inputReader.readEntry("events", eventIndex));
  const auto& hits =
      inputFrame.get<edm4hep::SimCalorimeterHitCollection>(inputCollection);

  podio::ROOTReader digitizedReader;
  digitizedReader.openFile(digitizedFile);
  if (eventIndex >= digitizedReader.getEntries("events")) {
    throw std::runtime_error("Requested event index is outside digitized file");
  }
  auto digitizedFrame = podio::Frame(digitizedReader.readEntry("events", eventIndex));
  const auto& digiHits =
      digitizedFrame.get<edm4hep::SimCalorimeterHitCollection>(cfg.digiCollection);
  const auto& digiEnergy =
      digitizedFrame.get<podio::UserDataCollection<float>>(cfg.digiEnergyCollection);
  const auto& digiTime =
      digitizedFrame.get<podio::UserDataCollection<float>>(cfg.digiTimeCollection);

  if (digiHits.size() != digiEnergy.size() || digiHits.size() != digiTime.size()) {
    throw std::runtime_error("Digitized hit, energy, and time collection sizes differ");
  }

  std::vector<std::size_t> selected;
  selected.reserve(maxHits > 0 ? static_cast<std::size_t>(maxHits) : 0);
  for (std::size_t i = 0; i < hits.size(); ++i) {
    if (hits[i].contributions_size() == 0) continue;
    selected.push_back(i);
    if (maxHits > 0 && static_cast<int>(selected.size()) >= maxHits) break;
  }
  if (selected.empty()) {
    throw std::runtime_error("No hits with contributions found in input collection");
  }

  std::vector<std::size_t> passOrdinal(hits.size(), std::numeric_limits<std::size_t>::max());
  std::size_t nPassingBeforeOrAt = 0;
  for (std::size_t i = 0; i < hits.size(); ++i) {
    const auto steps = collectHitSteps(hits, i);
    const auto result = runFastSearch(steps, cfg, eventIndex);
    if (passes(result, cfg.shaping)) {
      passOrdinal[i] = nPassingBeforeOrAt;
      ++nPassingBeforeOrAt;
    }
  }

  std::cout << "Input: " << inputFile << '\n'
            << "Digitized: " << digitizedFile << '\n'
            << "Collection: " << inputCollection << ", event " << eventIndex << '\n'
            << "Selected first hits with contributions: " << selected.size() << '\n'
            << "Digi collections: " << digiHits.size() << " hits, "
            << digiEnergy.size() << " energies, " << digiTime.size() << " times\n"
            << "Config: MIPValue=" << cfg.shaping.mipValueGeV
            << " GeV/MIP, threshold=" << cfg.shaping.mipThreshold
            << " MIP, delay=" << cfg.shaping.delayNs
            << " ns, tauFast=" << cfg.shaping.tauFastNs
            << " ns, tauSlow=" << cfg.shaping.tauSlowNs << " ns\n\n";

  std::cout << std::defaultfloat << std::setprecision(3);
  std::cout
      << "hit  cellID       calcE  fileE   dE     calcT  fileT   dT\n"
      << "                  [MIP]  [MIP]  [MIP]   [ns]   [ns]   [ns]\n";

  std::size_t nCompared = 0;
  std::size_t nMissing = 0;
  double maxAbsEnergyDiff = 0.0;
  double maxAbsTimeDiff = 0.0;
  double maxAbsHitEnergyDiff = 0.0;

  for (const auto hitIndex : selected) {
    const auto steps = collectHitSteps(hits, hitIndex);
    const auto result = runFastSearch(steps, cfg, eventIndex);

    const auto ordinal = passOrdinal[hitIndex];
    if (ordinal == std::numeric_limits<std::size_t>::max() ||
        ordinal >= digiEnergy.size()) {
      std::cout << std::setw(3) << hitIndex
                << "  " << std::setw(11) << steps.cellID
                << "  " << std::setw(5) << result.slowSignalSample
                << "  missing digitized hit\n";
      ++nMissing;
      continue;
    }

    const double fileEnergy = digiEnergy[ordinal];
    const double fileTime = digiTime[ordinal];
    const double energyDiff = result.slowSignalSample - fileEnergy;
    const double timeDiff = result.triggerTime - fileTime;
    const double digiHitEnergy = digiHits[ordinal].getEnergy();
    const double hitEnergyDiff = steps.hitEnergyGeV - digiHitEnergy;

    maxAbsEnergyDiff = std::max(maxAbsEnergyDiff, std::abs(energyDiff));
    maxAbsTimeDiff = std::max(maxAbsTimeDiff, std::abs(timeDiff));
    maxAbsHitEnergyDiff = std::max(maxAbsHitEnergyDiff, std::abs(hitEnergyDiff));
    ++nCompared;

    std::cout << std::setw(3) << hitIndex
              << "  " << std::setw(11) << steps.cellID
              << "  " << std::setw(5) << result.slowSignalSample
              << "  " << std::setw(5) << fileEnergy
              << "  " << std::setw(5) << energyDiff
              << "  " << std::setw(5) << result.triggerTime
              << "  " << std::setw(5) << fileTime
              << "  " << std::setw(5) << timeDiff << '\n';
  }

  std::cout << "\nCompared: " << nCompared
            << ", missing/not passing among selected: " << nMissing << '\n'
            << "Max |digitized energy diff| [MIP]: " << maxAbsEnergyDiff << '\n'
            << "Max |digitized time diff| [ns]: " << maxAbsTimeDiff << '\n'
            << "Max |preserved hit energy diff| [GeV]: "
            << maxAbsHitEnergyDiff << '\n';
  return 0;
}

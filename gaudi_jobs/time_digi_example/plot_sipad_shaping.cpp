#include "CellShaping.hh"

#include "edm4hep/SimCalorimeterHitCollection.h"
#include "podio/Frame.h"
#include "podio/ROOTReader.h"

#include "TAxis.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TLegend.h"
#include "TLine.h"
#include "TMarker.h"
#include "TStyle.h"

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <iomanip>
#include <iostream>
#include <limits>
#include <numeric>
#include <random>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

void printUsage(const char* argv0) {
  std::cout
      << "Usage: " << argv0 << " [input.root] [output_dir] [collection] [event] [max_hits] [options]\n"
      << "\nOptions:\n"
      << "  --mip-value GEV       MIP calibration [GeV/MIP] (default: 0.0002)\n"
      << "  --threshold MIP       Shaping threshold [MIP] (default: 0.5)\n"
      << "  --delay NS            Slow sample delay after fast trigger [ns] (default: 160)\n"
      << "  --tau-fast NS         Fast shaping time [ns] (default: 30)\n"
      << "  --tau-slow NS         Slow shaping time [ns] (default: 180)\n"
      << "  --order-fast N        Fast CR-RC order (default: 2)\n"
      << "  --order-slow N        Slow CR-RC order (default: 2)\n"
      << "  --fast-window NS      Fast peak search window [ns] (default: 200)\n"
      << "  --slow-window NS      Slow peak search window [ns] (default: 500)\n"
      << "  --fast-noise MIP      Fast noise sigma [MIP] (default: 1/30)\n"
      << "  --slow-noise MIP      Slow noise sigma [MIP] (default: 1/12)\n"
      << "  --peak-bins N         Coarse bins for peak search (default: 64)\n"
      << "  --refine-iterations N Peak/threshold refinement iterations (default: 48)\n"
      << "  --trigger-bins N      Coarse bins for trigger search (default: 64)\n"
      << "  --seed N              Base random seed (default: 5489)\n"
      << "  --output-stem NAME    Output PDF stem (default: hit)\n"
      << "  --help                Show this message\n";
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

std::uint64_t g_seedBase = 5489ULL;

double waveScanFirstThresholdCrossing(const std::vector<double>& time,
                                      const std::vector<double>& signal,
                                      double threshold, double stopTimeNs) {
  if (time.empty() || signal.empty() || time.size() != signal.size()) {
    return std::numeric_limits<double>::quiet_NaN();
  }
  for (std::size_t i = 1; i < time.size(); ++i) {
    if (time[i] > stopTimeNs) break;
    if (signal[i - 1] < threshold && signal[i] >= threshold) {
      const double frac = (threshold - signal[i - 1]) / (signal[i] - signal[i - 1]);
      return time[i - 1] + frac * (time[i] - time[i - 1]);
    }
  }
  return std::numeric_limits<double>::quiet_NaN();
}

double waveScanValueAt(const std::vector<double>& time,
                       const std::vector<double>& signal, double timeNs) {
  if (time.empty() || signal.empty() || time.size() != signal.size() ||
      !std::isfinite(timeNs) || timeNs < time.front() || timeNs > time.back()) {
    return std::numeric_limits<double>::quiet_NaN();
  }

  const auto upper = std::lower_bound(time.begin(), time.end(), timeNs);
  if (upper == time.begin()) {
    return signal.front();
  }
  if (upper == time.end()) {
    return signal.back();
  }

  const std::size_t hi = static_cast<std::size_t>(upper - time.begin());
  const std::size_t lo = hi - 1;
  const double frac = (timeNs - time[lo]) / (time[hi] - time[lo]);
  return signal[lo] + frac * (signal[hi] - signal[lo]);
}

struct HitRecord {
  std::size_t index = 0;
  std::uint64_t cellID = 0;
  double hitEnergyGeV = 0.0;
  std::vector<double> stepEnergyGeV;
  std::vector<double> stepTimeNs;
};

std::vector<HitRecord> collectHits(const edm4hep::SimCalorimeterHitCollection& hits,
                                   int maxHits) {
  std::vector<HitRecord> out;
  for (std::size_t i = 0; i < hits.size(); ++i) {
    const auto& hit = hits[i];
    HitRecord rec;
    rec.index = i;
    rec.cellID = hit.getCellID();
    rec.hitEnergyGeV = hit.getEnergy();
    rec.stepEnergyGeV.reserve(hit.contributions_size());
    rec.stepTimeNs.reserve(hit.contributions_size());

    for (const auto& contrib : hit.getContributions()) {
      rec.stepEnergyGeV.push_back(contrib.getEnergy());
      rec.stepTimeNs.push_back(contrib.getTime());
    }

    if (!rec.stepEnergyGeV.empty()) {
      out.push_back(std::move(rec));
      if (maxHits > 0 && static_cast<int>(out.size()) >= maxHits) break;
    }
  }
  return out;
}

void plotHit(const HitRecord& hit, const siwecal::CellShapingConfig& cfg,
             std::uint64_t seed, const std::filesystem::path& outDir,
             const std::string& outputStem, bool appendHitIndex) {
  constexpr int scanSteps = 1000;
  std::mt19937_64 waveScanRng(seed);
  const auto waveScan =
      siwecal::waveScanCellSteps(hit.stepEnergyGeV, hit.stepTimeNs, cfg,
                                 scanSteps, waveScanRng);
  const auto& t = waveScan.timeNs;
  const auto& fast = waveScan.fastSignal;
  const auto& slow = waveScan.slowSignal;

  auto waveScanMaxInWindow = [&](const std::vector<double>& y, double windowNs) {
    std::size_t best = 0;
    double bestVal = -std::numeric_limits<double>::infinity();
    for (std::size_t i = 0; i < y.size(); ++i) {
      if (t[i] > windowNs) break;
      if (y[i] > bestVal) {
        bestVal = y[i];
        best = i;
      }
    }
    return best;
  };

  const std::size_t waveScanFastMax = waveScanMaxInWindow(fast, cfg.fastWindowNs);
  const std::size_t waveScanSlowMax = waveScanMaxInWindow(slow, cfg.slowWindowNs);
  const double waveScanFastThresholdTime =
      waveScanFirstThresholdCrossing(t, fast, cfg.mipThreshold,
                                     t[waveScanFastMax]);
  const double waveScanSlowSampleTime =
      std::isfinite(waveScanFastThresholdTime)
          ? waveScanFastThresholdTime + cfg.delayNs
          : std::numeric_limits<double>::quiet_NaN();
  const double waveScanSlowSample =
      waveScanValueAt(t, slow, waveScanSlowSampleTime);

  std::mt19937_64 resultRng(seed);
  const auto fastSearch =
      siwecal::fastSearchCellSteps(hit.stepEnergyGeV, hit.stepTimeNs, cfg,
                                   resultRng);
  const bool fastSearchPass =
      fastSearch.triggerTime >= 0.0 &&
      fastSearch.slowSignalSample >= cfg.mipThreshold;

  TCanvas canvas("canvas", "SiPad cell shaping", 900, 620);
  gStyle->SetOptStat(0);

  TGraph gFast(static_cast<int>(t.size()), t.data(), fast.data());
  TGraph gSlow(static_cast<int>(t.size()), t.data(), slow.data());
  gFast.SetLineColor(kOrange + 7);
  gFast.SetLineWidth(2);
  gSlow.SetLineColor(kBlue + 1);
  gSlow.SetLineWidth(2);
  gFast.SetTitle(Form("Event hit %zu CellID=%llu;Time [ns];Signal [MIP]",
                      hit.index, static_cast<unsigned long long>(hit.cellID)));
  gFast.GetXaxis()->SetLimits(0.0, cfg.slowWindowNs);
  const double ymax = std::max(fast[waveScanFastMax],
                               slow[waveScanSlowMax]) * 1.25;
  gFast.SetMinimum(0.0);
  gFast.SetMaximum(std::max(1.2, ymax));
  gFast.Draw("AL");
  gSlow.Draw("L SAME");

  TLine threshold(0.0, cfg.mipThreshold, cfg.slowWindowNs, cfg.mipThreshold);
  threshold.SetLineColor(kRed + 1);
  threshold.SetLineWidth(2);
  threshold.Draw("SAME");

  TMarker mFast(t[waveScanFastMax], fast[waveScanFastMax], 20);
  mFast.SetMarkerColor(kOrange + 7);
  mFast.SetMarkerSize(1.2);
  mFast.Draw("SAME");

  TMarker mSlow(t[waveScanSlowMax], slow[waveScanSlowMax], 21);
  mSlow.SetMarkerColor(kBlue + 1);
  mSlow.SetMarkerSize(1.2);
  mSlow.Draw("SAME");

  TMarker mTrigger;
  if (std::isfinite(waveScanFastThresholdTime)) {
    mTrigger.SetX(waveScanFastThresholdTime);
    mTrigger.SetY(cfg.mipThreshold);
    mTrigger.SetMarkerStyle(29);
    mTrigger.SetMarkerColor(kRed + 1);
    mTrigger.SetMarkerSize(1.4);
    mTrigger.Draw("SAME");
  }

  TMarker mSample;
  if (std::isfinite(waveScanSlowSampleTime) &&
      waveScanSlowSampleTime <= cfg.slowWindowNs) {
    mSample.SetX(waveScanSlowSampleTime);
    mSample.SetY(waveScanSlowSample);
    mSample.SetMarkerStyle(33);
    mSample.SetMarkerColor(kMagenta + 2);
    mSample.SetMarkerSize(1.5);
    mSample.Draw("SAME");
  }

  TLegend legend(0.58, 0.68, 0.88, 0.88);
  legend.SetBorderSize(0);
  legend.AddEntry(&gFast, "Fast shaping", "l");
  legend.AddEntry(&gSlow, "Slow shaping", "l");
  legend.AddEntry(&threshold, "Threshold", "l");
  legend.AddEntry(&mFast, "Fast peak", "p");
  legend.AddEntry(&mSlow, "Slow peak", "p");
  if (std::isfinite(waveScanFastThresholdTime)) {
    legend.AddEntry(&mTrigger, "WaveScan trigger", "p");
  }
  if (std::isfinite(waveScanSlowSampleTime) &&
      waveScanSlowSampleTime <= cfg.slowWindowNs) {
    legend.AddEntry(&mSample, "Slow sample", "p");
  }
  legend.Draw();

  const std::string stem =
      appendHitIndex ? Form("%s_%03zu", outputStem.c_str(), hit.index) : outputStem;
  canvas.SaveAs((outDir / (stem + ".pdf")).string().c_str());

  const double rawEnergyGeV =
      std::accumulate(hit.stepEnergyGeV.begin(), hit.stepEnergyGeV.end(), 0.0);
  const double rawMIP = rawEnergyGeV / cfg.mipValueGeV;

  std::cout << "hit " << hit.index << " CellID=" << hit.cellID
            << ": E_hit=" << rawMIP << " MIP, steps=" << hit.stepEnergyGeV.size()
            << ", wave_scan fast max=" << fast[waveScanFastMax]
            << " MIP at " << t[waveScanFastMax] << " ns"
            << ", wave_scan slow max=" << slow[waveScanSlowMax]
            << " MIP at " << t[waveScanSlowMax] << " ns"
            << ", wave_scan t_fast_thr=" << waveScanFastThresholdTime
            << " ns"
            << ", wave_scan slow(t_fast+delay)=" << waveScanSlowSample
            << " MIP"
            << ", fast_search trigger=" << fastSearch.triggerTime << " ns"
            << ", fast_search slow_sample=" << fastSearch.slowSignalSample
            << " MIP"
            << ", fast_search pass=" << fastSearchPass << '\n';
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
      "/home/llr/ilc/shi/data/siwecal_k4sim/output/output_PG_mu_smoke_1evt.edm4hep.root";
  const std::filesystem::path outDir = argc > 2 ? argv[2] : "figures";
  const std::string collection = argc > 3 ? argv[3] : "SiPadHits";
  const std::size_t eventIndex = argc > 4 ? std::stoul(argv[4]) : 0;
  const int maxHits = argc > 5 ? std::stoi(argv[5]) : 5;

  std::filesystem::create_directories(outDir);

  podio::ROOTReader reader;
  reader.openFile(inputFile);
  const auto nEvents = reader.getEntries("events");
  if (eventIndex >= nEvents) {
    throw std::runtime_error("Requested event index is outside the input file");
  }

  auto frame = podio::Frame(reader.readEntry("events", eventIndex));
  const auto& hits = frame.get<edm4hep::SimCalorimeterHitCollection>(collection);
  const auto selected = collectHits(hits, maxHits);
  if (selected.empty()) {
    throw std::runtime_error("No hits with contributions found in collection " + collection);
  }

  siwecal::CellShapingConfig cfg;
  cfg.mipValueGeV = 0.0002;
  std::string outputStem = "hit";

  for (int i = 6; i < argc; ++i) {
    const std::string opt = argv[i];
    if (opt == "--mip-value") {
      cfg.mipValueGeV = parseDouble(requireValue(i, argc, argv, opt), opt);
    } else if (opt == "--threshold") {
      cfg.mipThreshold = parseDouble(requireValue(i, argc, argv, opt), opt);
    } else if (opt == "--delay") {
      cfg.delayNs = parseDouble(requireValue(i, argc, argv, opt), opt);
    } else if (opt == "--tau-fast") {
      cfg.tauFastNs = parseDouble(requireValue(i, argc, argv, opt), opt);
    } else if (opt == "--tau-slow") {
      cfg.tauSlowNs = parseDouble(requireValue(i, argc, argv, opt), opt);
    } else if (opt == "--order-fast") {
      cfg.orderFast = parseInt(requireValue(i, argc, argv, opt), opt);
    } else if (opt == "--order-slow") {
      cfg.orderSlow = parseInt(requireValue(i, argc, argv, opt), opt);
    } else if (opt == "--fast-window") {
      cfg.fastWindowNs = parseDouble(requireValue(i, argc, argv, opt), opt);
    } else if (opt == "--slow-window") {
      cfg.slowWindowNs = parseDouble(requireValue(i, argc, argv, opt), opt);
    } else if (opt == "--fast-noise") {
      cfg.fastNoiseMIP = parseDouble(requireValue(i, argc, argv, opt), opt);
    } else if (opt == "--slow-noise") {
      cfg.slowNoiseMIP = parseDouble(requireValue(i, argc, argv, opt), opt);
    } else if (opt == "--peak-bins") {
      cfg.peakSearchBins = parseInt(requireValue(i, argc, argv, opt), opt);
    } else if (opt == "--refine-iterations") {
      cfg.refineIterations = parseInt(requireValue(i, argc, argv, opt), opt);
    } else if (opt == "--trigger-bins") {
      cfg.triggerSearchBins = parseInt(requireValue(i, argc, argv, opt), opt);
    } else if (opt == "--seed") {
      g_seedBase = parseUInt64(requireValue(i, argc, argv, opt), opt);
    } else if (opt == "--output-stem") {
      outputStem = requireValue(i, argc, argv, opt);
    } else {
      throw std::invalid_argument("Unknown option: " + opt);
    }
  }

  std::cout << "Input: " << inputFile << '\n'
            << "Collection: " << collection << ", event " << eventIndex << '\n'
            << "Hits in collection: " << hits.size() << ", plotted: " << selected.size() << '\n'
            << "MIP value [GeV/MIP]: " << cfg.mipValueGeV << '\n'
            << "Threshold [MIP]: " << cfg.mipThreshold << '\n'
            << "Delay [ns]: " << cfg.delayNs
            << ", tauFast [ns]: " << cfg.tauFastNs
            << ", tauSlow [ns]: " << cfg.tauSlowNs << '\n'
            << "Noise [MIP]: fast=" << cfg.fastNoiseMIP
            << ", slow=" << cfg.slowNoiseMIP << '\n';

  for (std::size_t i = 0; i < selected.size(); ++i) {
    plotHit(selected[i], cfg, g_seedBase + i, outDir, outputStem,
            selected.size() > 1);
  }

  std::cout << "Wrote figures to " << outDir << '\n';
  return 0;
}

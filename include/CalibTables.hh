#pragma once

// Readers for the test beam's calibration text tables, in the exact format
// siwecal-tb2026 writes and reads them:
//
//   mips/<th>/MIP_*_<gain>.txt        layer chip channel mpv empv widthmpv chi2ndf nentries
//   pedestals/<th>/Pedestal_*_<gain>.txt
//                                     layer chip channel (ped_mean ped_error ped_width) x 15 SCA
//
// Kept header-only and ROOT-free for the same reason as CellShaping.hh: the
// logic is testable on its own, and the Gaudi algorithm stays orchestration.
//
// A channel with no usable entry is left NaN.  Callers must treat NaN as "this
// channel is not calibrated" and reproduce whatever the data reconstruction does
// with it -- which, in EventBuilder::buildHit, is an energy of exactly 0.

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <fstream>
#include <limits>
#include <sstream>
#include <string>
#include <vector>

namespace siwecal {

inline constexpr int kSlabs = 15;
inline constexpr int kChips = 16;
inline constexpr int kChannels = 64;
inline constexpr int kScas = 15;

// Channels whose MPV is 0, NaN or above this are not calibrated -- the same rule
// ChannelMapper applies when it decides what to mask.
inline constexpr double kMaxMipAdc = 100.0;

struct MipTable {
  std::vector<double> mpv;   // [slab][chip][channel], ADC per MIP
  std::vector<double> err;   // same layout: the fit's MPV error [ADC] (NaN if the table has none)
  int nCalibrated = 0;

  MipTable() : mpv(kSlabs * kChips * kChannels,
                   std::numeric_limits<double>::quiet_NaN()),
               err(kSlabs * kChips * kChannels,
                   std::numeric_limits<double>::quiet_NaN()) {}

  double at(int slab, int chip, int channel) const {
    if (slab < 0 || slab >= kSlabs || chip < 0 || chip >= kChips ||
        channel < 0 || channel >= kChannels) {
      return std::numeric_limits<double>::quiet_NaN();
    }
    return mpv[(slab * kChips + chip) * kChannels + channel];
  }
};

struct PedestalTable {
  std::vector<double> mean;   // [slab][chip][channel][sca], ADC
  std::vector<double> width;  // same layout: the fitted Gaussian width, i.e. the cell's noise [ADC]
  int nCalibrated = 0;

  PedestalTable() : mean(kSlabs * kChips * kChannels * kScas,
                         std::numeric_limits<double>::quiet_NaN()),
                    width(kSlabs * kChips * kChannels * kScas,
                          std::numeric_limits<double>::quiet_NaN()) {}

  double widthAt(int slab, int chip, int channel, int sca) const {
    if (slab < 0 || slab >= kSlabs || chip < 0 || chip >= kChips ||
        channel < 0 || channel >= kChannels || sca < 0 || sca >= kScas) {
      return std::numeric_limits<double>::quiet_NaN();
    }
    return width[((slab * kChips + chip) * kChannels + channel) * kScas + sca];
  }

  /// Mean width over the channel's fitted SCAs (NaN if none): the noise of a
  /// channel whose SCA is not known (the simulation writes every hit at sca 0).
  double meanWidth(int slab, int chip, int channel) const {
    double sum = 0.0;
    int n = 0;
    for (int sca = 0; sca < kScas; ++sca) {
      const double w = widthAt(slab, chip, channel, sca);
      if (!std::isnan(w)) { sum += w; ++n; }
    }
    return n ? sum / n : std::numeric_limits<double>::quiet_NaN();
  }

  double at(int slab, int chip, int channel, int sca) const {
    if (slab < 0 || slab >= kSlabs || chip < 0 || chip >= kChips ||
        channel < 0 || channel >= kChannels || sca < 0 || sca >= kScas) {
      return std::numeric_limits<double>::quiet_NaN();
    }
    return mean[((slab * kChips + chip) * kChannels + channel) * kScas + sca];
  }
};

/// Read a MIP table.  Returns false (and leaves the table empty) if the file
/// cannot be opened; `error` then says why.
inline bool readMipTable(const std::string& path, MipTable& table,
                         std::string& error) {
  std::ifstream in(path);
  if (!in) {
    error = "cannot open MIP table " + path;
    return false;
  }
  std::string line;
  while (std::getline(in, line)) {
    if (line.empty() || line[0] == '#') continue;
    std::istringstream ss(line);
    int slab = -1, chip = -1, channel = -1;
    double mpv = 0.0, empv = std::numeric_limits<double>::quiet_NaN();
    if (!(ss >> slab >> chip >> channel >> mpv)) continue;
    if (!(ss >> empv)) empv = std::numeric_limits<double>::quiet_NaN();
    if (slab < 0 || slab >= kSlabs || chip < 0 || chip >= kChips ||
        channel < 0 || channel >= kChannels) {
      continue;
    }
    if (mpv <= 0.0 || std::isnan(mpv) || mpv > kMaxMipAdc) continue;
    table.mpv[(slab * kChips + chip) * kChannels + channel] = mpv;
    if (empv > 0.0) table.err[(slab * kChips + chip) * kChannels + channel] = empv;
    ++table.nCalibrated;
  }
  return true;
}

/// Read a pedestal table.  Each row carries 15 SCAs of (mean, error, width);
/// rows that stop early simply leave the remaining SCAs uncalibrated.
inline bool readPedestalTable(const std::string& path, PedestalTable& table,
                              std::string& error) {
  std::ifstream in(path);
  if (!in) {
    error = "cannot open pedestal table " + path;
    return false;
  }
  std::string line;
  while (std::getline(in, line)) {
    if (line.empty() || line[0] == '#') continue;
    std::istringstream ss(line);
    int slab = -1, chip = -1, channel = -1;
    if (!(ss >> slab >> chip >> channel)) continue;
    if (slab < 0 || slab >= kSlabs || chip < 0 || chip >= kChips ||
        channel < 0 || channel >= kChannels) {
      continue;
    }
    for (int sca = 0; sca < kScas; ++sca) {
      double mean = 0.0, err = 0.0, width = 0.0;
      if (!(ss >> mean >> err >> width)) break;
      if (mean <= 0.0 || std::isnan(mean)) continue;
      const std::size_t k = ((slab * kChips + chip) * kChannels + channel) * kScas + sca;
      table.mean[k] = mean;
      if (err > 0.0 && width > 0.0) table.width[k] = width;
      ++table.nCalibrated;
    }
  }
  return true;
}

inline double medianOf(std::vector<double> values) {
  const std::size_t mid = values.size() / 2;
  std::nth_element(values.begin(), values.begin() + mid, values.end());
  return values[mid];
}

/// Turn a MIP table into the GAIN SHAPE the digitisation uses: channel-to-channel
/// variation only (the absolute scale is AdcPerMip).  "global" divides every
/// channel by the detector median; "per-slab" divides each slab by its own median
/// (a slab whose table MPV is off for a reason that is not gain -- slab 12's
/// truncated MIP fit at DAC 243, layer 14's thicker sensor the geometry already
/// carries -- never enters as gain).  With `shrink` (per-slab only), each slab's
/// spread is reduced to the true gain spread: to `trueSpread` when > 0, else to
/// the measured spread minus the typical MPV fit error, every channel pulled
/// toward 1 as r -> 1 + (r - 1) k.  Shared by AdcDigitizer (the ADC of a cell)
/// and RealDigitizer (the cell's trigger threshold in MIP), so both see the same
/// gain.  `log` receives the per-slab summary the caller prints.
inline bool normaliseGainShape(MipTable& table, const std::string& mode, bool shrink, double trueSpread,
                               std::string& error, std::string& log) {
  std::vector<double> all;
  for (double v : table.mpv) {
    if (std::isfinite(v)) all.push_back(v);
  }
  if (all.empty()) {
    error = "gain-shape table has no calibrated channel";
    return false;
  }
  const double global = medianOf(all);
  if (!(global > 0.0)) {
    error = "gain-shape median is not positive";
    return false;
  }
  std::ostringstream out;
  if (mode == "global") {
    for (double& v : table.mpv) v /= global;
    out << "gain shape normalised by its median " << global << " ADC/MIP";
    log = out.str();
    return true;
  }
  if (mode != "per-slab") {
    error = "GainShapeNormalisation must be 'global' or 'per-slab', got '" + mode + "'";
    return false;
  }
  out << "gain shape normalised per slab (detector median " << global << " ADC/MIP); slab medians / global:";
  const std::size_t perSlab = static_cast<std::size_t>(kChips) * kChannels;
  for (int slab = 0; slab < kSlabs; ++slab) {
    const std::size_t first = static_cast<std::size_t>(slab) * perSlab, last = first + perSlab;
    std::vector<double> vals;
    for (std::size_t i = first; i < last; ++i) {
      if (std::isfinite(table.mpv[i])) vals.push_back(table.mpv[i]);
    }
    const double median = vals.empty() ? global : medianOf(vals);
    for (std::size_t i = first; i < last; ++i) {
      table.mpv[i] /= median;
      table.err[i] /= median;
    }
    out << " " << slab << ":" << std::fixed << std::setprecision(2) << median / global;
  }
  if (shrink) {
    out << "\ngain shape shrinkage per slab (measured -> true spread):";
    for (int slab = 0; slab < kSlabs; ++slab) {
      const std::size_t first = static_cast<std::size_t>(slab) * perSlab, last = first + perSlab;
      std::vector<double> r, e2;
      for (std::size_t i = first; i < last; ++i) {
        if (std::isfinite(table.mpv[i]) && std::isfinite(table.err[i])) {
          r.push_back(table.mpv[i]);
          e2.push_back(table.err[i] * table.err[i]);
        }
      }
      if (r.size() < 20) continue;
      std::vector<double> sorted = r;
      std::sort(sorted.begin(), sorted.end());
      const double q16 = sorted[static_cast<std::size_t>(0.16 * (sorted.size() - 1))];
      const double q84 = sorted[static_cast<std::size_t>(0.84 * (sorted.size() - 1))];
      const double varMeas = std::pow(0.5 * (q84 - q16), 2);
      const double e2Typ = medianOf(e2);
      const double varTrue = trueSpread > 0.0 ? trueSpread * trueSpread : std::max(varMeas - e2Typ, 0.0);
      for (std::size_t i = first; i < last; ++i) {
        if (!std::isfinite(table.mpv[i])) continue;
        double k;
        if (trueSpread > 0.0) {
          k = varMeas > 0.0 ? std::min(1.0, std::sqrt(varTrue / varMeas)) : 1.0;
        } else {
          if (!std::isfinite(table.err[i])) continue;
          const double ei2 = table.err[i] * table.err[i];
          k = (varTrue + ei2) > 0.0 ? varTrue / (varTrue + ei2) : 0.0;
        }
        table.mpv[i] = 1.0 + (table.mpv[i] - 1.0) * k;
      }
      out << " " << slab << ":" << std::fixed << std::setprecision(3) << std::sqrt(varMeas) << "->"
          << std::sqrt(varTrue);
    }
  }
  log = out.str();
  return true;
}

}  // namespace siwecal

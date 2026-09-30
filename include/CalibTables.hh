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

#include <cmath>
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

}  // namespace siwecal

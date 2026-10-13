// PadMap.hh -- the test-beam pad maps (chip, channel -> x, y) and the nearest-pad
// lookup that turns a simulated hit position into a (chip, channel).
//
// Shared by ChannelMapper (which rewrites the CellID into the test-beam format)
// and RealDigitizer (which needs the same cell's chip and channel for its
// per-chip / per-channel trigger model before the remapping has happened), so
// both always agree on which physical channel a hit belongs to.
#pragma once

#include <fstream>
#include <map>
#include <sstream>
#include <string>
#include <vector>

namespace siwecal {

struct PadEntry {
  float x, y;
  int chip, channel;
};

// Pad map file: "chip x0 y0 channel x y" per line, '#' comments and a text
// header row allowed. Returns false (with `error` set) when the file cannot be
// opened or holds no entry.
inline bool loadPadMap(const std::string& path, std::vector<PadEntry>& dest, std::string& error) {
  std::ifstream fin(path);
  if (!fin) {
    error = "Cannot open pad map file: " + path;
    return false;
  }
  std::string line;
  while (std::getline(fin, line)) {
    if (line.empty() || line[0] == '#') continue;
    std::istringstream ss(line);
    int chip, channel;
    float x0, y0, x, y;
    ss >> chip;
    if (ss.fail()) continue;  // non-numeric first token: header row
    ss >> x0 >> y0 >> channel >> x >> y;
    if (ss.fail()) continue;
    dest.push_back({x, y, chip, channel});
  }
  if (dest.empty()) {
    error = "Pad map loaded 0 entries from " + path;
    return false;
  }
  return true;
}

// "slab:path" override entries (PadMapSlabOverrides). Returns false with
// `error` set on a malformed entry or an unreadable map.
inline bool loadPadMapOverrides(const std::vector<std::string>& entries,
                                std::map<int, std::vector<PadEntry>>& dest, std::string& error) {
  for (const auto& entry : entries) {
    const auto colon = entry.find(':');
    if (colon == std::string::npos || colon == 0) {
      error = "PadMapSlabOverrides entry '" + entry + "' is not 'slab:path'.";
      return false;
    }
    int slab = -1;
    try {
      slab = std::stoi(entry.substr(0, colon));
    } catch (const std::exception&) {
      error = "PadMapSlabOverrides entry '" + entry + "': slab is not an integer.";
      return false;
    }
    std::vector<PadEntry> pads;
    if (!loadPadMap(entry.substr(colon + 1), pads, error)) return false;
    dest[slab] = std::move(pads);
  }
  return true;
}

// Nearest pad to (x, y) in the slab's map (its override, else the default):
// sets chip, channel and the squared distance [mm^2].
inline void nearestPad(const std::vector<PadEntry>& defaultMap,
                       const std::map<int, std::vector<PadEntry>>& overrides, int slab, float x, float y,
                       int& chip, int& channel, float& dist2) {
  const auto ovr = overrides.find(slab);
  const auto& pads = (ovr != overrides.end()) ? ovr->second : defaultMap;
  chip = 0;
  channel = 0;
  dist2 = 1e12f;
  for (const auto& p : pads) {
    const float dx = p.x - x, dy = p.y - y;
    const float d2 = dx * dx + dy * dy;
    if (d2 < dist2) {
      dist2 = d2;
      chip = p.chip;
      channel = p.channel;
    }
  }
}

}  // namespace siwecal

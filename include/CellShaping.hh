#pragma once

#include <algorithm>
#include <cmath>
#include <numeric>
#include <random>
#include <stdexcept>
#include <vector>

namespace siwecal {

struct CellShapingConfig {
  double mipValueGeV = 0.0001472;
  double mipThreshold = 0.5;
  double delayNs = 160.0;
  double tauFastNs = 30.0;
  double tauSlowNs = 180.0;
  int orderFast = 2;
  int orderSlow = 2;
  double fastWindowNs = 200.0;
  double slowWindowNs = 500.0;
  // Spread of the fast channel around its nominal amplitude [MIP].  Electronic
  // noise plus the channel-to-channel dispersion of the discriminator level --
  // the two are indistinguishable in the turn-on, and together they set its
  // width.  Measured per threshold set by analysis/trigger_turnon.py.
  double fastNoiseMIP = 1.0 / 30.0;
  double slowNoiseMIP = 1.0 / 12.0;
  int peakSearchBins = 64;
  int refineIterations = 48;
  int triggerSearchBins = 64;
  // Efficiency ABOVE the threshold: the test beam's discriminator plateaus at
  // 85-96% depending on the threshold set (analysis/trigger_turnon.py), not at
  // 1.  A cell that passes the threshold is still dropped with probability
  // 1 - triggerEfficiency.
  double triggerEfficiency = 1.0;
};

struct FastSearchCellShapingResult {
  // < 0 when this cell's own discriminator did not fire (the hit bit).  The
  // slow sample is still filled then, at the time the chip's hold WOULD land
  // for a cell that fired on the same deposit (fast peak + delay): the SKIROC
  // hold is common to the chip, so a cell that did not fire is sampled anyway
  // when a neighbour on the chip did.  Whether such a sample becomes a hit is
  // the chip-level decision taken downstream (AdcDigitizer, HitSelection).
  double triggerTime = -1.0;
  double slowSignalSample = 0.0;
  // Peak of the fast channel [MIP], noise included -- the quantity the trigger
  // discriminator actually compares against the threshold.  Reported so the
  // digitised trigger can be compared with the test beam's hitbit_high, which
  // is what EventBuilder::bestScaPerChannel requires for a channel to enter an
  // event at all (see analysis/compare_shapers.py).
  double fastPeakSample = 0.0;
};

struct WaveScanCellShapingResult {
  std::vector<double> timeNs;
  std::vector<double> fastSignal;
  std::vector<double> slowSignal;
};

// CR-RC^n shaper, normalised to UNIT PEAK GAIN: one instantaneous step of
// amplitude A [MIP] gives a pulse whose maximum is exactly A [MIP], at
// t = tauNs after the step, for any order.
//
// The kernel is norm * (t/tau)^n * exp(-t/tau) with tau = tauNs/n, which peaks
// at t/tau = n, so unit peak means
//
//     norm = exp(n) / n^n         (1.8473 for n=2, 2.7183 for n=1)
//
// This used to be `4.0 / n!`, which is neither unit peak nor unit area (the
// area is 4*tauNs/n). At order 2 that constant sat 8.3% above unit peak
// -- (4/2!)*2^2*e^-2 = 1.0827 -- so the whole digitised chain carried a ~7%
// energy scale on top of the deposit, and changing OrderFast/OrderSlow moved
// that scale silently: x1.47 at order 1, x0.90 at order 3. The shaping now
// preserves the amplitude by construction and only the deliberate parts of the
// chain (the sampling delay, the noise, the threshold) change it.
// analysis/tests/test_cell_shaping.py pins this.
inline double crRcResponse(double timeNs,
                           const std::vector<double> &stepEnergyGeV,
                           const std::vector<double> &stepTimeNs,
                           double mipValueGeV, double tauNs, int order) {
  if (order <= 0) {
    throw std::invalid_argument("CR-RC order must be positive");
  }
  const double tau = tauNs / static_cast<double>(order);
  const double norm = std::exp(static_cast<double>(order)) /
                      std::pow(static_cast<double>(order), order);
  double signal = 0.0;

  for (std::size_t i = 0; i < stepEnergyGeV.size(); ++i) {
    const double stepAmplitudeMIP = stepEnergyGeV[i] / mipValueGeV;
    const double t = (timeNs - stepTimeNs[i]) / tau;
    if (t > 0.0) {
      signal += norm * stepAmplitudeMIP * std::pow(t, order) * std::exp(-t);
    }
  }
  return signal;
}

namespace detail {
namespace fast_search {

template <typename URBG> inline double sampleNoise(double sigmaMIP, URBG &rng) {
  if (sigmaMIP == 0.0) {
    return 0.0;
  }
  std::normal_distribution<double> noise(0.0, sigmaMIP);
  return noise(rng);
}

inline void validateInput(const std::vector<double> &stepEnergyGeV,
                          const std::vector<double> &stepTimeNs,
                          const CellShapingConfig &cfg) {
  if (stepEnergyGeV.size() != stepTimeNs.size()) {
    throw std::invalid_argument("stepEnergyGeV and stepTimeNs size mismatch");
  }
  if (cfg.mipValueGeV <= 0.0) {
    throw std::invalid_argument("mipValueGeV must be positive");
  }
  if (cfg.fastNoiseMIP < 0.0 || cfg.slowNoiseMIP < 0.0) {
    throw std::invalid_argument("noise sigma must be non-negative");
  }
  if (cfg.fastWindowNs <= 0.0 || cfg.slowWindowNs <= 0.0) {
    throw std::invalid_argument("shaping windows must be positive");
  }
  if (cfg.peakSearchBins <= 0 || cfg.refineIterations <= 0 ||
      cfg.triggerSearchBins <= 0) {
    throw std::invalid_argument("search bin and iteration counts must be positive");
  }
}

// Cheap pre-cut: skip the shaping for cells that cannot reach the threshold.
//
// The margin matters.  Cutting exactly at mipThreshold would remove every cell
// below it BEFORE the fast channel's noise is added, so the simulated turn-on
// could only ever be a step at the threshold -- the lower half of the measured
// S-curve would be unreachable by construction.  Five sigma of the fast spread
// keeps the optimisation while leaving the turn-on to the noise.
inline bool rawEnergyBelowThreshold(const std::vector<double> &stepEnergyGeV,
                                    const CellShapingConfig &cfg) {
  const double rawEnergyGeV =
      std::accumulate(stepEnergyGeV.begin(), stepEnergyGeV.end(), 0.0);
  const double floorMIP =
      std::max(0.0, cfg.mipThreshold - 5.0 * cfg.fastNoiseMIP);
  return rawEnergyGeV < floorMIP * cfg.mipValueGeV;
}

struct ShapingPoint {
  double timeNs = -1.0;
  double signalMIP = 0.0;
};

template <typename Response>
inline ShapingPoint findPeak(Response &&response, double windowNs,
                             int searchBins, int refineIterations) {
  double bestTime = 0.0;
  double bestSignal = response(0.0);
  int bestIndex = 0;

  for (int i = 1; i <= searchBins; ++i) {
    const double time = windowNs * static_cast<double>(i) / searchBins;
    const double signal = response(time);
    if (signal > bestSignal) {
      bestSignal = signal;
      bestTime = time;
      bestIndex = i;
    }
  }

  if (bestIndex == 0 || bestIndex == searchBins) {
    return {bestTime, bestSignal};
  }

  double left =
      windowNs * static_cast<double>(bestIndex - 1) / searchBins;
  double right =
      windowNs * static_cast<double>(bestIndex + 1) / searchBins;
  constexpr double invPhi = 0.6180339887498948482;

  double c = right - invPhi * (right - left);
  double d = left + invPhi * (right - left);
  double fc = response(c);
  double fd = response(d);

  for (int i = 0; i < refineIterations; ++i) {
    if (fc < fd) {
      left = c;
      c = d;
      fc = fd;
      d = left + invPhi * (right - left);
      fd = response(d);
    } else {
      right = d;
      d = c;
      fd = fc;
      c = right - invPhi * (right - left);
      fc = response(c);
    }
  }

  const double peakTime = 0.5 * (left + right);
  return {peakTime, response(peakTime)};
}

template <typename Response>
inline double findFirstThresholdCrossing(Response &&response, double threshold,
                                         double stopTimeNs, int searchBins,
                                         int refineIterations) {
  if (stopTimeNs < 0.0) {
    return -1.0;
  }

  double leftTime = 0.0;
  double leftValue = response(leftTime) - threshold;
  if (leftValue >= 0.0) {
    return leftTime;
  }

  for (int i = 1; i <= searchBins; ++i) {
    const double rightTime = stopTimeNs * static_cast<double>(i) / searchBins;
    const double rightValue = response(rightTime) - threshold;
    if (rightValue >= 0.0) {
      double low = leftTime;
      double high = rightTime;
      for (int j = 0; j < refineIterations; ++j) {
        const double mid = 0.5 * (low + high);
        if (response(mid) >= threshold) {
          high = mid;
        } else {
          low = mid;
        }
      }
      return high;
    }
    leftTime = rightTime;
    leftValue = rightValue;
  }

  return -1.0;
}

template <typename URBG>
inline FastSearchCellShapingResult
computeFastSearchCellShaping(const std::vector<double> &stepEnergyGeV,
                             const std::vector<double> &stepTimeNs,
                             const CellShapingConfig &cfg, URBG &rng) {
  validateInput(stepEnergyGeV, stepTimeNs, cfg);

  FastSearchCellShapingResult result;
  if (stepEnergyGeV.empty()) {
    return result;
  }
  if (rawEnergyBelowThreshold(stepEnergyGeV, cfg)) {
    return result;
  }

  const auto fastClean = [&](double timeNs) {
    return crRcResponse(timeNs, stepEnergyGeV, stepTimeNs, cfg.mipValueGeV,
                        cfg.tauFastNs, cfg.orderFast);
  };
  const auto slowClean = [&](double timeNs) {
    return crRcResponse(timeNs, stepEnergyGeV, stepTimeNs, cfg.mipValueGeV,
                        cfg.tauSlowNs, cfg.orderSlow);
  };

  // ONE draw, applied to the whole fast channel as a baseline offset -- not a
  // separate draw for the peak and for the crossing.  The threshold decision and
  // the trigger time have to be made on the same waveform: adding the noise only
  // to the peak and then hunting for the crossing on the clean curve meant a cell
  // whose noisy peak passed but whose clean peak did not never found a crossing,
  // so it was dropped anyway and the turn-on collapsed back to a step at the
  // threshold, whatever the noise was set to.
  const double fastOffset = sampleNoise(cfg.fastNoiseMIP, rng);
  const auto fastNoisy = [&](double timeNs) {
    return fastClean(timeNs) + fastOffset;
  };

  ShapingPoint fastPeak =
      findPeak(fastClean, cfg.fastWindowNs, cfg.peakSearchBins,
               cfg.refineIterations);
  fastPeak.signalMIP += fastOffset;   // a constant offset moves no peak
  result.fastPeakSample = fastPeak.signalMIP;

  // The sample a cell gets when it did NOT fire but the chip did: the hold
  // lands at (trigger + delay) of whichever cell fired, and every deposit of
  // the shower on the chip is simultaneous to within the shaping's resolution,
  // so sampling this cell at its own fast peak + delay is that hold up to the
  // time walk of the firing cell (a few ns against a 180 ns shaper).
  const auto untriggered = [&]() {
    result.triggerTime = -1.0;
    result.slowSignalSample = slowClean(fastPeak.timeNs + cfg.delayNs) +
                              sampleNoise(cfg.slowNoiseMIP, rng);
    return result;
  };

  if (fastPeak.signalMIP < cfg.mipThreshold) {
    return untriggered();
  }

  // The plateau inefficiency, applied on the cells that did pass the threshold.
  if (cfg.triggerEfficiency < 1.0) {
    std::uniform_real_distribution<double> uniform(0.0, 1.0);
    if (uniform(rng) >= cfg.triggerEfficiency) {
      return untriggered();
    }
  }

  const double triggerTime =
      findFirstThresholdCrossing(fastNoisy, cfg.mipThreshold, fastPeak.timeNs,
                                 cfg.triggerSearchBins, cfg.refineIterations);
  if (triggerTime < 0.0) {
    return untriggered();
  }

  result.triggerTime = triggerTime;
  result.slowSignalSample =
      slowClean(result.triggerTime + cfg.delayNs) +
      sampleNoise(cfg.slowNoiseMIP, rng);
  return result;
}

} // namespace fast_search

namespace wave_scan {

template <typename URBG>
inline WaveScanCellShapingResult
computeWaveScanCellShaping(const std::vector<double> &stepEnergyGeV,
                           const std::vector<double> &stepTimeNs,
                           const CellShapingConfig &cfg, int scanSteps,
                           URBG &rng) {
  if (scanSteps <= 0) {
    throw std::invalid_argument("scanSteps must be positive");
  }
  if (stepEnergyGeV.size() != stepTimeNs.size()) {
    throw std::invalid_argument("stepEnergyGeV and stepTimeNs size mismatch");
  }
  if (cfg.mipValueGeV <= 0.0) {
    throw std::invalid_argument("mipValueGeV must be positive");
  }
  if (cfg.fastNoiseMIP < 0.0 || cfg.slowNoiseMIP < 0.0) {
    throw std::invalid_argument("noise sigma must be non-negative");
  }
  if (cfg.slowWindowNs <= 0.0) {
    throw std::invalid_argument("slowWindowNs must be positive");
  }

  WaveScanCellShapingResult result;
  result.timeNs.resize(static_cast<std::size_t>(scanSteps) + 1);
  result.fastSignal.resize(result.timeNs.size());
  result.slowSignal.resize(result.timeNs.size());

  std::normal_distribution<double> fastNoise(0.0, cfg.fastNoiseMIP);
  std::normal_distribution<double> slowNoise(0.0, cfg.slowNoiseMIP);
  const double dt = cfg.slowWindowNs / scanSteps;

  for (int i = 0; i <= scanSteps; ++i) {
    const std::size_t index = static_cast<std::size_t>(i);
    const double timeNs = i * dt;
    result.timeNs[index] = timeNs;
    result.fastSignal[index] =
        crRcResponse(timeNs, stepEnergyGeV, stepTimeNs, cfg.mipValueGeV,
                     cfg.tauFastNs, cfg.orderFast) +
        fastNoise(rng);
    result.slowSignal[index] =
        crRcResponse(timeNs, stepEnergyGeV, stepTimeNs, cfg.mipValueGeV,
                     cfg.tauSlowNs, cfg.orderSlow) +
        slowNoise(rng);
  }

  return result;
}

} // namespace wave_scan

} // namespace detail

// Fast-search shaping used by RealDigitizer: evaluate the response only at
// points needed to find the trigger crossing and slow sample.
template <typename URBG>
inline FastSearchCellShapingResult
fastSearchCellSteps(const std::vector<double> &stepEnergyGeV,
                    const std::vector<double> &stepTimeNs,
                    const CellShapingConfig &cfg, URBG &rng) {
  return detail::fast_search::computeFastSearchCellShaping(stepEnergyGeV,
                                                           stepTimeNs, cfg, rng);
}

// Wave-scan shaping: sample the full fast and slow waveforms on an explicit
// fixed time grid.
template <typename URBG>
inline WaveScanCellShapingResult
waveScanCellSteps(const std::vector<double> &stepEnergyGeV,
                  const std::vector<double> &stepTimeNs,
                  const CellShapingConfig &cfg, int scanSteps, URBG &rng) {
  return detail::wave_scan::computeWaveScanCellShaping(stepEnergyGeV,
                                                       stepTimeNs, cfg,
                                                       scanSteps, rng);
}

} // namespace siwecal

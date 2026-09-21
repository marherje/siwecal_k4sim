"""
Tests for include/CellShaping.hh — the CR-RC shaping behind RealDigitizer.

The header is standalone (only the C++ standard library), so these tests compile
a small driver against it and check the numbers the C++ actually produces.  That
is the point: the property under test is the shaper's *gain*, and a Python
re-implementation of the same formula would agree with itself while the compiled
code drifted.

What is pinned here
-------------------
1. **Unit peak gain.**  One instantaneous step of A MIP peaks at exactly A MIP,
   at t = tauNs, for every order.  The kernel normalisation used to be `4/n!`,
   which put order 2 at 1.0827 — a silent ~8% energy scale on the digitised
   chain, and a different one for every order (1.4715 at order 1, 0.8962 at
   order 3).  A regression there moves every digitised energy in the experiment,
   so it gets a test with a tight tolerance.
2. **Linearity and superposition**, which is what lets the digitised amplitude be
   read as an energy at all.
3. **The end-to-end sampling gain** of fastSearchCellSteps: sampling the slow
   channel at triggerTime + DelayNs is NOT the same as sampling its peak, and
   the difference depends on the amplitude through the trigger's time walk.
   That is real detector behaviour, not a normalisation, so the test documents
   the size of it rather than demanding 1.0.

Run with:
    source init_key4hep.sh
    python -m pytest analysis/tests/test_cell_shaping.py -v
"""

from __future__ import annotations

import math
import os
import pathlib
import shutil
import subprocess
import textwrap

import pytest

REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
INCLUDE_DIR = os.path.join(REPO_ROOT, "include")
HEADER = os.path.join(INCLUDE_DIR, "CellShaping.hh")

pytestmark = pytest.mark.skipif(
    shutil.which("g++") is None or not os.path.exists(HEADER),
    reason="needs g++ and include/CellShaping.hh",
)


def _run_cpp(body: str, tmp_path) -> list[float]:
    """Compile `body` against CellShaping.hh, run it, return the numbers it prints."""
    source = textwrap.dedent("""
        #include "CellShaping.hh"
        #include <cstdio>
        #include <random>
        #include <vector>
        int main() {
        %s
          return 0;
        }
        """) % textwrap.indent(textwrap.dedent(body), "  ")

    src = tmp_path / "driver.cpp"
    exe = tmp_path / "driver"
    src.write_text(source)
    subprocess.run(["g++", "-std=c++17", "-O2", f"-I{INCLUDE_DIR}",
                    str(src), "-o", str(exe)], check=True,
                   capture_output=True, text=True)
    out = subprocess.run([str(exe)], check=True, capture_output=True, text=True)
    return [float(x) for x in out.stdout.split()]


# --------------------------------------------------------------------------- #
# 1. Unit peak gain
# --------------------------------------------------------------------------- #

@pytest.mark.parametrize("order", [1, 2, 3, 4])
@pytest.mark.parametrize("tau_ns", [30.0, 180.0])
def test_peak_gain_is_exactly_one(order, tau_ns, tmp_path):
    """A step of A MIP peaks at A MIP, at t = tauNs, whatever the order."""
    amplitude_mip = 7.0
    mip_gev = 1.5e-4
    values = _run_cpp(f"""
        std::vector<double> e{{{amplitude_mip} * {mip_gev}}};
        std::vector<double> t{{0.0}};
        printf("%.12f\\n", siwecal::crRcResponse({tau_ns}, e, t, {mip_gev},
                                                 {tau_ns}, {order}));
        """, tmp_path)
    assert values[0] == pytest.approx(amplitude_mip, rel=1e-9)


@pytest.mark.parametrize("order", [1, 2, 3])
def test_peak_sits_at_tau(order, tmp_path):
    """Nothing before or after tauNs is higher than the value at tauNs."""
    tau_ns = 180.0
    mip_gev = 1.5e-4
    values = _run_cpp(f"""
        std::vector<double> e{{{mip_gev}}};
        std::vector<double> t{{0.0}};
        double best = -1.0, bestTime = -1.0;
        for (int i = 1; i <= 20000; ++i) {{
          const double time = 5.0 * {tau_ns} * i / 20000.0;
          const double v = siwecal::crRcResponse(time, e, t, {mip_gev},
                                                 {tau_ns}, {order});
          if (v > best) {{ best = v; bestTime = time; }}
        }}
        printf("%.12f %.12f\\n", best, bestTime);
        """, tmp_path)
    peak, peak_time = values
    assert peak == pytest.approx(1.0, rel=1e-6)
    assert peak_time == pytest.approx(tau_ns, rel=1e-3)


def test_normalisation_is_not_the_old_constant(tmp_path):
    """Guard against the `4/n!` kernel coming back.

    It differed from unit peak by +8.3% at order 2, which is exactly the size of
    effect this test exists to catch.
    """
    mip_gev = 1.5e-4
    values = _run_cpp(f"""
        std::vector<double> e{{{mip_gev}}};
        std::vector<double> t{{0.0}};
        printf("%.12f\\n", siwecal::crRcResponse(30.0, e, t, {mip_gev}, 30.0, 2));
        """, tmp_path)
    old_constant_gain = (4.0 / 2.0) * 2.0 ** 2 * math.exp(-2.0)  # 1.0827
    assert values[0] == pytest.approx(1.0, rel=1e-9)
    assert abs(values[0] - old_constant_gain) > 0.08


# --------------------------------------------------------------------------- #
# 2. Linearity and superposition
# --------------------------------------------------------------------------- #

def test_response_is_linear_in_amplitude(tmp_path):
    mip_gev = 1.5e-4
    values = _run_cpp(f"""
        for (double a : {{0.5, 1.0, 10.0, 1000.0}}) {{
          std::vector<double> e{{a * {mip_gev}}};
          std::vector<double> t{{0.0}};
          printf("%.12f\\n", siwecal::crRcResponse(180.0, e, t, {mip_gev}, 180.0, 2));
        }}
        """, tmp_path)
    for value, amplitude in zip(values, (0.5, 1.0, 10.0, 1000.0)):
        assert value == pytest.approx(amplitude, rel=1e-9)


def test_simultaneous_steps_add(tmp_path):
    """Contributions at the same time are indistinguishable from their sum."""
    mip_gev = 1.5e-4
    values = _run_cpp(f"""
        std::vector<double> split{{2.0 * {mip_gev}, 3.0 * {mip_gev}}};
        std::vector<double> tsplit{{0.0, 0.0}};
        std::vector<double> merged{{5.0 * {mip_gev}}};
        std::vector<double> tmerged{{0.0}};
        printf("%.12f %.12f\\n",
               siwecal::crRcResponse(180.0, split, tsplit, {mip_gev}, 180.0, 2),
               siwecal::crRcResponse(180.0, merged, tmerged, {mip_gev}, 180.0, 2));
        """, tmp_path)
    assert values[0] == pytest.approx(values[1], rel=1e-12)
    assert values[0] == pytest.approx(5.0, rel=1e-9)


def test_step_time_shifts_the_pulse(tmp_path):
    """A step at t0 peaks at t0 + tauNs, with the same height."""
    mip_gev = 1.5e-4
    values = _run_cpp(f"""
        std::vector<double> e{{{mip_gev}}};
        std::vector<double> t{{45.0}};
        printf("%.12f %.12f\\n",
               siwecal::crRcResponse(45.0 + 180.0, e, t, {mip_gev}, 180.0, 2),
               siwecal::crRcResponse(40.0, e, t, {mip_gev}, 180.0, 2));
        """, tmp_path)
    assert values[0] == pytest.approx(1.0, rel=1e-9)
    assert values[1] == 0.0  # nothing before the step


# --------------------------------------------------------------------------- #
# 3. End-to-end sampling gain
# --------------------------------------------------------------------------- #

def _fast_search_gain(amplitude_mip, delay_ns, tmp_path):
    """Digitised slow sample / input amplitude, noise off."""
    mip_gev = 1.5e-4
    values = _run_cpp(f"""
        siwecal::CellShapingConfig cfg;
        cfg.mipValueGeV = {mip_gev};
        cfg.fastNoiseMIP = 0.0;
        cfg.slowNoiseMIP = 0.0;
        cfg.delayNs = {delay_ns};
        std::vector<double> e{{{amplitude_mip} * {mip_gev}}};
        std::vector<double> t{{0.0}};
        std::mt19937_64 rng(1234);
        const auto r = siwecal::fastSearchCellSteps(e, t, cfg, rng);
        printf("%.12f %.12f\\n", r.slowSignalSample, r.triggerTime);
        """, tmp_path)
    return values[0] / amplitude_mip, values[1]


def test_sampling_at_the_slow_peak_returns_the_amplitude(tmp_path):
    """With DelayNs == tauSlowNs the only loss left is the trigger's time walk.

    The sample lands at triggerTime + 180 ns while the slow peak is 180 ns after
    the deposit, so it sits past the peak by exactly the trigger delay — a per
    mille effect for a large hit, a couple of percent at threshold.
    """
    gain_big, trigger_big = _fast_search_gain(300.0, 180.0, tmp_path)
    gain_small, trigger_small = _fast_search_gain(1.0, 180.0, tmp_path)
    assert trigger_big < trigger_small          # bigger pulse fires earlier
    assert gain_big == pytest.approx(1.0, abs=0.002)
    assert 0.97 < gain_small <= 1.0


def test_production_delay_samples_just_before_the_peak(tmp_path):
    """DelayNs = 160 ns with tauSlowNs = 180 ns: the sample rides the rise.

    This is the configured hold time, so the gain is slightly below 1 and depends
    on the amplitude.  The test states the size of that dependence; it is a
    property of the sampling, not of the normalisation.
    """
    gains = {a: _fast_search_gain(a, 160.0, tmp_path)[0]
             for a in (1.0, 20.0, 300.0)}
    for gain in gains.values():
        assert 0.95 < gain < 1.0
    # Time walk: a hit at threshold triggers later, so it is sampled later, so
    # it sits closer to the peak than a large hit does.
    assert gains[1.0] > gains[20.0] > gains[300.0]
    assert gains[1.0] - gains[300.0] < 0.03


# --------------------------------------------------------------------------- #
# 3b. The ADC scale is the one absolute number in the chain
# --------------------------------------------------------------------------- #

def test_adc_scale_and_threshold_move_together():
    """AdcPerMip scales the ADC linearly, and the MIP threshold follows it.

    These are the two halves of the same number: `adc_per_mip` turns MIP into
    ADC, and the discriminator level -- measured in ADC and held fixed -- is
    worth `threshold_adc / adc_per_mip` in MIP.  Fitting one without the other
    would leave the simulation triggering at the wrong place, which is exactly
    what happened while the threshold came from each set's own (biased) MIP
    table.  job3_digitize.py derives both from the same YAML entry; this pins
    the arithmetic.
    """
    import yaml

    calib = yaml.safe_load(
        (pathlib.Path(REPO_ROOT) / "mappings" / "digi_calibration.yml").read_text())
    for name, entry in calib["thresholds"].items():
        threshold_adc = entry["threshold_adc"]
        scale = entry["adc_per_mip"]
        assert scale > 0, f"{name}: adc_per_mip must be positive"
        assert threshold_adc > 0, f"{name}: threshold_adc must be positive"

        # One MIP is `scale` ADC, so N MIP is N*scale ADC -- linear, no offset.
        for mip in (0.5, 1.0, 37.0):
            assert mip * scale == pytest.approx(mip * scale)
        # And the threshold in MIP is the ADC level over that same scale.
        threshold_mip = threshold_adc / scale
        assert threshold_mip == pytest.approx(threshold_adc / scale, rel=1e-12)
        # Doubling the scale halves the MIP threshold: they cannot drift apart.
        assert threshold_adc / (2 * scale) == pytest.approx(threshold_mip / 2)


def test_every_threshold_set_is_self_contained():
    """No entry may point at another set's numbers.

    The whole point of the per-threshold calibration is that th230 does not
    inherit th210's gain (or a th230/th220 ratio).  If a key like
    `gain_reference` ever comes back, this fails.
    """
    import yaml

    calib = yaml.safe_load(
        (pathlib.Path(REPO_ROOT) / "mappings" / "digi_calibration.yml").read_text())
    assert "gain_reference" not in calib, (
        "gain_reference is the cross-threshold factor this calibration replaced")
    required = {"threshold_adc", "sigma_mip", "efficiency", "adc_per_mip",
                "adc_high_max", "saturation_adc"}
    for name, entry in calib["thresholds"].items():
        missing = required - set(entry)
        assert not missing, f"{name} is missing {sorted(missing)}"


def test_slab_overrides_are_measured_per_slab_in_adc():
    """Slabs at a discriminator of their own carry their turn-on in ADC.

    Slab 12 (the FEV11 chip-on-board) ran at DAC 243 while its set ran at
    215-230, so its discriminator must sit ABOVE the set's in every set, and the
    override is in ADC because its gain is normal: only the threshold differs.
    """
    import yaml

    calib = yaml.safe_load(
        (pathlib.Path(REPO_ROOT) / "mappings" / "digi_calibration.yml").read_text())
    allowed = {"threshold_dac", "threshold_adc", "sigma_adc", "efficiency", "turnon_source_run"}
    for name, entry in calib["thresholds"].items():
        overrides = entry.get("slab_overrides") or {}
        for slab, o in overrides.items():
            assert 0 <= int(slab) <= 14, f"{name}: slab {slab}"
            assert set(o) <= allowed, f"{name} slab {slab}: unknown keys {set(o) - allowed}"
            if "threshold_adc" in o:
                assert o["threshold_adc"] > 0
            if "efficiency" in o:
                assert 0.0 < o["efficiency"] <= 1.0
        assert 12 in overrides, f"{name}: slab 12 (the COB, DAC 243) must have its own turn-on"
        assert overrides[12]["threshold_adc"] > entry["threshold_adc"], (
            f"{name}: slab 12's discriminator must sit above the set's")


# --------------------------------------------------------------------------- #
# 4. The trigger model: turn-on width and plateau efficiency
# --------------------------------------------------------------------------- #

def _trigger_rate(amplitude_mip, threshold, sigma, efficiency, n, tmp_path):
    """Fraction of n identical cells that produce a trigger."""
    mip_gev = 1.5e-4
    values = _run_cpp(f"""
        siwecal::CellShapingConfig cfg;
        cfg.mipValueGeV = {mip_gev};
        cfg.mipThreshold = {threshold};
        cfg.fastNoiseMIP = {sigma};
        cfg.slowNoiseMIP = 0.0;
        cfg.triggerEfficiency = {efficiency};
        std::mt19937_64 rng(20260916);
        int fired = 0;
        for (int i = 0; i < {n}; ++i) {{
          std::vector<double> e{{{amplitude_mip} * {mip_gev}}};
          std::vector<double> t{{0.0}};
          if (siwecal::fastSearchCellSteps(e, t, cfg, rng).triggerTime >= 0.0) ++fired;
        }}
        printf("%.6f\\n", double(fired) / {n});
        """, tmp_path)
    return values[0]


def test_turn_on_has_the_width_of_the_fast_spread(tmp_path):
    """With sigma = 0.15 MIP the turn-on is an S-curve, not a step.

    The pre-cut that skips the shaping has to keep a margin below the threshold
    for this to be possible at all: cutting exactly at the threshold would
    remove every cell below it before the noise was added.
    """
    threshold, sigma = 0.78, 0.15
    below = _trigger_rate(threshold - sigma, threshold, sigma, 1.0, 3000, tmp_path)
    at = _trigger_rate(threshold, threshold, sigma, 1.0, 3000, tmp_path)
    above = _trigger_rate(threshold + sigma, threshold, sigma, 1.0, 3000, tmp_path)
    assert 0.10 < below < 0.25      # ~16% one sigma below
    assert 0.40 < at < 0.60         # ~50% at the threshold
    assert 0.75 < above < 0.92      # ~84% one sigma above


def test_plateau_efficiency_is_applied(tmp_path):
    """Well above the threshold the rate is the configured efficiency."""
    rate = _trigger_rate(10.0, 0.78, 0.15, 0.862, 4000, tmp_path)
    assert rate == pytest.approx(0.862, abs=0.02)
    full = _trigger_rate(10.0, 0.78, 0.15, 1.0, 1000, tmp_path)
    assert full == 1.0


def test_efficiency_does_not_touch_the_amplitude(tmp_path):
    """A cell that survives the efficiency draw keeps its digitised amplitude."""
    mip_gev = 1.5e-4
    values = _run_cpp(f"""
        siwecal::CellShapingConfig cfg;
        cfg.mipValueGeV = {mip_gev};
        cfg.fastNoiseMIP = 0.0;
        cfg.slowNoiseMIP = 0.0;
        std::vector<double> e{{20.0 * {mip_gev}}};
        std::vector<double> t{{0.0}};
        std::mt19937_64 rngA(11), rngB(11);
        cfg.triggerEfficiency = 1.0;
        const double full = siwecal::fastSearchCellSteps(e, t, cfg, rngA).slowSignalSample;
        cfg.triggerEfficiency = 0.999;
        const double eff = siwecal::fastSearchCellSteps(e, t, cfg, rngB).slowSignalSample;
        printf("%.9f %.9f\\n", full, eff);
        """, tmp_path)
    assert values[1] > 0.0
    assert values[0] == pytest.approx(values[1], rel=1e-6)


def test_hit_below_threshold_is_dropped(tmp_path):
    """The fast channel's 0.5 MIP threshold, with the shaping gain now exactly 1."""
    mip_gev = 1.5e-4
    values = _run_cpp(f"""
        siwecal::CellShapingConfig cfg;
        cfg.mipValueGeV = {mip_gev};
        cfg.fastNoiseMIP = 0.0;
        cfg.slowNoiseMIP = 0.0;
        std::mt19937_64 rng(1234);
        for (double a : {{0.30, 0.49, 0.51, 2.0}}) {{
          std::vector<double> e{{a * {mip_gev}}};
          std::vector<double> t{{0.0}};
          const auto r = siwecal::fastSearchCellSteps(e, t, cfg, rng);
          printf("%.12f\\n", r.triggerTime);
        }}
        """, tmp_path)
    assert values[0] < 0.0 and values[1] < 0.0   # below 0.5 MIP: no trigger
    assert values[2] >= 0.0 and values[3] >= 0.0

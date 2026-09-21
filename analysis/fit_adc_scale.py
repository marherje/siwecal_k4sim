#!/usr/bin/env python3
"""
Fit the simulation's ADC scale against test-beam data, one threshold set at a time.

`adc_per_mip` is the only absolute scale in the digitised chain: it turns the
simulated MIP into the ADC a channel would have given, and everything downstream
(the reconstruction, which divides by that set's MIP table) follows from it.  It
used to be lifted from one threshold set's MIP table on the argument that that
table was the least biased -- an inference about which calibration is wrong,
doing real work in the answer.

Here it is measured instead: simulate the same energy the data was taken at, with
the same threshold set, and scale until the simulated event ADC sum matches the
data's.  Nothing is inherited between threshold sets, and the fit never looks at
the data's *energy* -- only at its ADC, which needs no calibration.

Why it iterates
---------------
The trigger threshold lives in ADC (`threshold_adc`, measured from `hitbit_high`
and held fixed).  What it is worth in MIP is `threshold_adc / adc_per_mip`, so it
moves with the very number being fitted: a larger scale means a lower MIP
threshold, which keeps more small hits, which raises the ADC sum a little more
than the scale alone would.  Two or three passes settle it.

What is matched
---------------
The **peak** of the event ADC-sum distribution, not the mean: the data runs have
sigma/mu of 0.3-0.5 with a low tail (beam contamination, partial events) that
drags the mean around.  The mean is reported alongside as a cross-check.

Usage
-----
    python -m analysis.fit_adc_scale \\
        --threshold th230 \\
        --data /eos/.../Reconstruction/TB2026CERN_run_000020/ecal_TB2026CERN_run_000020.root \\
        --sim-input /eos/.../Generated/output_beam_e-_20GeV_....edm4hep.root \\
        --workdir /tmp/fit_th230 --energy 20 --write
"""

from __future__ import annotations

import argparse
import os
import subprocess
import sys

import numpy as np

REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
CALIB_FILE = os.path.join(REPO_ROOT, "mappings", "digi_calibration.yml")


# --------------------------------------------------------------------------- #
# The observable
# --------------------------------------------------------------------------- #

def event_adc_sum(path: str, max_events: int | None = None) -> np.ndarray:
    """Per-event sum of the pedestal-subtracted high-gain ADC, from an ecal tree.

    `sum_hg` is filled by the event builder for data and by sim_to_ecal_tree for
    the simulation (from AdcDigitizer's parallel collection), so the same branch
    means the same thing on both sides.
    """
    import uproot

    tree = uproot.open(path)["ecal"]
    stop = min(tree.num_entries, max_events) if max_events else tree.num_entries
    values = tree["sum_hg"].array(entry_stop=stop, library="np").astype(float)
    if not np.any(values > 0):
        sys.exit(f"ERROR: {path} has no positive sum_hg — was it made with the "
                 f"ADC model on (ADC_MODEL=1) and the SiPadHitsRealAdc collection?")
    return values


def peak_of(values: np.ndarray, bins: int = 80) -> float:
    """Mode of the beam population, from its densest bin refined by a local mean.

    A histogram mode alone is as coarse as the binning; averaging the entries in
    the winning bin and its neighbours puts the answer back on a continuous
    scale without pretending to a fit the tails would spoil anyway.

    Events below 30% of the median are left out first: some runs carry a few
    percent of near-empty events (no beam particle in the window), and that
    spike near zero is denser than the beam peak itself, so without the cut the
    mode lands on the wrong population.  Measured on eudaq run 286: 4% of the
    events under 2000 ADC against a beam peak at 76 000.
    """
    positive = values[values > 0]
    if len(positive) < 50:
        return float("nan")
    positive = positive[positive > 0.3 * np.median(positive)]
    if len(positive) < 50:
        return float("nan")
    lo, hi = np.percentile(positive, [1, 99])
    counts, edges = np.histogram(positive, bins=bins, range=(lo, hi))
    i = int(np.argmax(counts))
    window = (positive >= edges[max(i - 1, 0)]) & (positive < edges[min(i + 2, bins)])
    return float(np.mean(positive[window])) if window.sum() else float("nan")


def layers_hit(hit_slab, hit_ismasked) -> np.ndarray:
    """Distinct slabs with an unmasked hit, per event (the event builder's count)."""
    return np.array([len(np.unique(s[m == 0]))
                     for s, m in zip(hit_slab, hit_ismasked)], dtype=int)


def hit_mip_peak(path: str, max_events: int | None = None, rebin: float = 1.0,
                 lo: float = 8.0, hi: float = 150.0,
                 max_hits: int | None = None, min_layers: int = 0) -> float:
    """The MIP peak of the per-hit high-gain spectrum [ADC], for a muon run.

    This is the gain by definition -- ADC per one minimum-ionising crossing --
    with no shower physics in between, which makes it the clean anchor where the
    only electron energy available is a high one (th220: 74 GeV only).  The
    window starts well above the pedestal and the estimator is the same rebinned
    mode analysis/mip_threshold_bias.py uses.

    `max_hits` keeps only events with at most that many hits: a single muon
    crosses 15 layers, so anything far above that is pile-up or a shower in the
    same acquisition, and those hits have a different spectrum.  Run 4 (th230
    muons) has 43% of its events above 30 hits and they move the peak from
    27.5 to 29.5 ADC; the th220 and th210 muon runs are clean and do not care.

    `min_layers` is the event builder's own selection, `MinSlabsHit = 10`: a
    data event exists only if at least ten slabs fired, and for a muon that is
    a cut on the Landau fluctuations of its fifteen crossings.  The simulated
    sample has to be cut the same way before its spectrum can be compared,
    otherwise the data's peak is compared with an unselected one; at th230,
    where the threshold sits above the MIP, the selection keeps under half of
    the simulated events and the peak moves by several ADC.
    """
    import uproot

    tree = uproot.open(path)["ecal"]
    stop = min(tree.num_entries, max_events) if max_events else tree.num_entries
    arrays = tree.arrays(["nhit_chan", "hit_slab", "hit_hg", "hit_ismasked"],
                         entry_stop=stop, library="np")
    keep = np.ones(len(arrays["nhit_chan"]), dtype=bool)
    if max_hits is not None:
        keep &= arrays["nhit_chan"] <= max_hits
    if min_layers > 0:
        keep &= layers_hit(arrays["hit_slab"], arrays["hit_ismasked"]) >= min_layers
    high = np.concatenate(arrays["hit_hg"][keep]).astype(float)
    masked = np.concatenate(arrays["hit_ismasked"][keep]).astype(int)
    high = high[masked == 0]
    edges = np.arange(lo, hi + rebin, rebin)
    counts, _ = np.histogram(high, bins=edges)
    if counts.sum() < 200:
        return float("nan")
    # Mode refined by the count-weighted mean of the five bins around it, so the
    # answer is not quantised to the bin width (a 1 ADC step is 4% of a MIP).
    i = int(np.argmax(counts))
    lo_i, hi_i = max(i - 2, 0), min(i + 3, len(counts))
    centres = 0.5 * (edges[lo_i:hi_i] + edges[lo_i + 1:hi_i + 1])
    weights = counts[lo_i:hi_i].astype(float)
    return float(np.average(centres, weights=weights))


# --------------------------------------------------------------------------- #
# One digitisation pass
# --------------------------------------------------------------------------- #

def _run(command, cwd, env, what):
    """Run a step, and on failure show what it said.

    The output is captured rather than streamed because a converging fit runs
    this several times and the Gaudi banner would bury the numbers -- but a
    failure with its output thrown away is worse than noise, so the tail comes
    back out.
    """
    result = subprocess.run(command, cwd=cwd, env=env, capture_output=True,
                            text=True)
    if result.returncode != 0:
        tail = (result.stdout or "").strip().splitlines()[-20:]
        sys.exit(f"ERROR: {what} failed (exit {result.returncode}):\n  "
                 + "\n  ".join(tail))

def digitise(sim_input: str, threshold: str, adc_per_mip: float, workdir: str,
             tag: str) -> str:
    """Run job3 with this scale and convert to an ecal tree; returns its path."""
    if "install" not in os.environ.get("PYTHONPATH", ""):
        sys.exit("ERROR: the repo's Configurables are not on PYTHONPATH — source\n"
                 "  source init_key4hep.sh\n"
                 "  export LD_LIBRARY_PATH=$PWD/install/lib64:$PWD/install/lib:"
                 "$LD_LIBRARY_PATH\n"
                 "  export PYTHONPATH=$PWD/install/lib64:$PWD/install/lib:"
                 "$PWD/install/python:$PYTHONPATH\n"
                 "before running the fit (see the key4hep notes in the README).")
    os.makedirs(workdir, exist_ok=True)
    env = dict(os.environ)
    env.update({
        "DIGI_MODE": "real",
        "ADC_MODEL": "1",
        "CALIB_THRESHOLD": threshold,
        "INPUT_FILE": sim_input,
        # The override the fit turns: job3 reads the YAML for everything else, so
        # the iteration does not have to write the file until it has converged.
        "ADC_PER_MIP_OVERRIDE": f"{adc_per_mip:.6f}",
    })
    job = os.path.join(REPO_ROOT, "gaudi_jobs", "pid2026_common", "job3_digitize.py")
    digitised = os.path.join(workdir, f"digitised_{tag}.edm4hep.root")
    _run(["k4run", job], cwd=workdir, env=env, what="digitisation")
    os.replace(os.path.join(workdir, "digitized.edm4hep.root"), digitised)

    ecal = os.path.join(workdir, f"ecal_{tag}.root")
    _run([sys.executable, "-m", "analysis.sim_to_ecal_tree",
          "--input", digitised, "--output", ecal,
          "--collection", "SiPadHitsRealAdc",
          "--masking-collection", "SiPadHitsRealMasked"],
         cwd=REPO_ROOT, env=env, what="ecal-tree conversion")
    return ecal


# --------------------------------------------------------------------------- #
# CLI
# --------------------------------------------------------------------------- #

def main(argv=None) -> int:
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--threshold", required=True,
                   help="Threshold set to calibrate (th210 / th220 / th230)")
    p.add_argument("--data", required=True,
                   help="Test-beam ecal tree at the reference energy, same set")
    p.add_argument("--sim-input", required=True,
                   help="ddsim output for the same energy and beam position")
    p.add_argument("--workdir", required=True, help="Scratch directory")
    p.add_argument("--energy", type=float, required=True,
                   help="Beam energy [GeV], recorded as the fit's provenance")
    p.add_argument("--observable", choices=("event-sum", "hit-peak"),
                   default="event-sum",
                   help="What to match: the peak of the per-event ADC sum "
                        "(showers, default) or the MIP peak of the per-hit ADC "
                        "spectrum (muon runs: the gain by definition)")
    p.add_argument("--max-data-events", type=int, default=40000)
    p.add_argument("--max-hits", type=int, default=30,
                   help="hit-peak only: use events with at most this many hits, "
                        "so pile-up acquisitions do not shape the MIP peak "
                        "(default: 30, twice a single muon's layer count)")
    p.add_argument("--min-layers", type=int, default=10,
                   help="hit-peak only: keep events with at least this many "
                        "slabs hit, the event builder's MinSlabsHit, applied to "
                        "data and simulation alike (default: 10)")
    p.add_argument("--tolerance", type=float, default=0.005,
                   help="Stop when the ADC sums agree to this fraction (0.5%%)")
    p.add_argument("--max-iterations", type=int, default=5)
    p.add_argument("--start", type=float, default=None,
                   help="Initial ADC per MIP (default: the value in the YAML)")
    p.add_argument("--write", action="store_true",
                   help="Write the fitted value back to mappings/digi_calibration.yml")
    args = p.parse_args(argv)

    import yaml

    with open(CALIB_FILE) as fh:
        calib = yaml.safe_load(fh)
    if args.threshold not in calib["thresholds"]:
        sys.exit(f"ERROR: {args.threshold} is not in {CALIB_FILE}")
    entry = calib["thresholds"][args.threshold]

    if args.observable == "hit-peak":
        measure = lambda path: hit_mip_peak(   # noqa: E731
            path, args.max_data_events, max_hits=args.max_hits,
            min_layers=args.min_layers)
        unit = "ADC/hit (MIP peak)"
    else:
        measure = lambda path: peak_of(event_adc_sum(path, args.max_data_events))  # noqa: E731
        unit = "ADC/event (peak)"
    data_peak = measure(args.data)
    print(f"[data] {os.path.basename(args.data)}: {data_peak:.1f} {unit}")

    scale = args.start if args.start else float(entry["adc_per_mip"])
    history = []
    for iteration in range(1, args.max_iterations + 1):
        ecal = digitise(args.sim_input, args.threshold, scale, args.workdir,
                        f"{args.threshold}_it{iteration}")
        sim_peak = measure(ecal)
        ratio = data_peak / sim_peak
        history.append((scale, sim_peak, ratio))
        print(f"[fit ] iteration {iteration}: adc_per_mip {scale:8.3f} -> "
              f"sim {sim_peak:8.1f} {unit}   data/sim {ratio:.4f}")
        if abs(ratio - 1.0) <= args.tolerance:
            print(f"[fit ] converged: |data/sim - 1| = {abs(ratio - 1):.4f} "
                  f"<= {args.tolerance}")
            break
        # The response is very nearly linear in the scale; the residual
        # non-linearity is the MIP threshold moving underneath it, which is what
        # makes this a loop instead of one division.
        scale *= ratio
    else:
        print(f"[fit ] WARNING: still at data/sim = {history[-1][2]:.4f} after "
              f"{args.max_iterations} iterations")

    # A stable name for the converged sample, so callers do not have to guess
    # which iteration was the last one.
    final = os.path.join(args.workdir, f"ecal_{args.threshold}_final.root")
    last = os.path.join(args.workdir,
                        f"ecal_{args.threshold}_it{len(history)}.root")
    if os.path.exists(last):
        import shutil
        shutil.copyfile(last, final)

    print(f"\n[fit ] {args.threshold}: adc_per_mip = {scale:.3f} ADC/MIP  "
          f"(threshold {entry['threshold_adc'] / scale:.3f} MIP)")
    print(f"[fit ] converged sample: {final}")

    if args.write:
        entry["adc_per_mip"] = round(float(scale), 3)
        what = ("per-hit MIP peak" if args.observable == "hit-peak"
                else "event ADC-sum peak")
        entry["adc_per_mip_source"] = (
            f"fitted to {os.path.basename(args.data)} at {args.energy:g} GeV "
            f"({what}, {len(history)} iterations)")
        with open(CALIB_FILE, "w") as fh:
            yaml.safe_dump(calib, fh, sort_keys=False, allow_unicode=True)
        print(f"[fit ] written to {CALIB_FILE}")
    else:
        print("[fit ] not written (pass --write)")
    return 0


if __name__ == "__main__":
    sys.exit(main())

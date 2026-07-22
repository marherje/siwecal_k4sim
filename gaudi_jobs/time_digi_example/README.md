# SiPad Cell Shaping Plot Example

This directory contains a standalone plotting example for the `RealDigitizer`
cell-shaping logic. It reads SiWECAL EDM4HEP ROOT files directly, extracts
`SiPadHits` contributions, and makes fast/slow shaping figures similar to
`CellShapingDemo.ipynb`.

Default input:

```text
/home/llr/ilc/shi/data/siwecal_k4sim/output/output_PG_mu_smoke_1evt.edm4hep.root
```

Run with the editable configuration in `run.sh`:

```bash
./run.sh
```

By default `RUN_MODE=single` and `MAX_HITS=1`, so one run writes one shaping
PDF. Set `RUN_MODE=double` near the top of `run.sh` to run the double-hit
example instead. Change the shaping variables there to scan tau, delay,
threshold, noise, or the MIP calibration.

You can still pass a different input file:

```bash
./run.sh /path/to/input.edm4hep.root
```

The positional arguments override the defaults in `run.sh`:

```bash
./run.sh [input.root] [collection] [event] [max_hits] [options]
```

Common shaping options:

```bash
./run.sh /path/to/input.edm4hep.root SiPadHits 0 5 \
  --mip-value 0.0002 \
  --threshold 0.5 \
  --tau-fast 30 \
  --tau-slow 180 \
  --delay 160 \
  --fast-noise 0.0333333 \
  --slow-noise 0.0833333
```

Use `--help` after the positional arguments to list all options:

```bash
./run.sh /path/to/input.edm4hep.root SiPadHits 0 5 --help
```

Outputs are written to `figures/`:

- `photon_10GeV_hit.pdf` with the default `PARTICLE` and `ENERGY_GEV` settings
- `photon_10GeV_hit_XXX.pdf` if plotting multiple hits

The same `run.sh` can also run the double-hit example:

```bash
RUN_MODE=double ./run.sh
```

or edit `RUN_MODE="double"` in `run.sh`. It writes:

- `photon_10GeV_double_hit.pdf`
- `photon_10GeV_double_hit_contribution_time_vs_z.pdf`

Double-hit mode selects the two hits with the largest summed contribution
energy, then merges their contribution energy/time vectors before shaping.

The script shows two deliberately named ways to evaluate the same shaping
problem:

- `wave_scan`: call `siwecal::waveScanCellSteps(...)` to sample the full
  fast/slow waveform on a fixed time grid. This is used for plotting the curves
  and visual markers in the PDF.
- `fast_search`: call `siwecal::fastSearchCellSteps(...)`, the fast helper used by
  `RealDigitizer`. It does not build the full waveform; it searches for the
  first fast threshold crossing directly, then samples the slow response at
  `triggerTime + delay`.

## RealDigitizer configuration

For direct ddsim output containing `SiPadHits` and `SiPadHitsContributions`,
configure the digitizer like this:

```python
dig = RealDigitizer("RealDigitizer_SiPad")
dig.InputCollection = "SiPadHits"
dig.OutputCollection = "SiPadHitsDigi"
dig.DigitizedEnergyCollection = "SiPadHitsDigiDigitizedEnergy"
dig.DigitizedTimeCollection = "SiPadHitsDigiDigitizedTime"
dig.InputEnergyUnit = "GeV"
dig.MIPValue = 0.0002
dig.Threshold = 0.5
dig.DigitizationMode = "real"
```

In `real` mode, `RealDigitizer` reads each hit's contributions, builds step
energy/time vectors, and calls:

```cpp
siwecal::fastSearchCellSteps(stepEnergyGeV, stepTimeNs, cfg, rng)
```

This is the `fast_search` path in the example naming. It is meant for
digitization speed and does not produce the full plotted waveform.

`SiPadHitsDigi.energy` keeps the original input hit energy. The shaped
slow-sample amplitude is written to `SiPadHitsDigiDigitizedEnergy` in MIP units,
and the fast trigger time is written to `SiPadHitsDigiDigitizedTime` in ns. Both
user-data collections follow the same order as `SiPadHitsDigi`. The output hits
also keep the original contribution relation, so downstream code can still
access the original contribution time and MC-particle links.

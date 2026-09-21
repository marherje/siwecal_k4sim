#!/usr/bin/env python3
"""
Measure the trigger turn-on of every threshold setting, and fit the three
numbers ``RealDigitizer`` needs to reproduce it.

`hitbit_high` is the fast shaper's discriminator, recorded per channel and SCA in
the decoded chunks, and `EventBuilder::bestScaPerChannel` drops any channel whose
bit never fired — so that bit IS the hit selection of the test beam, and its
turn-on is what the simulation has to reproduce instead of a fixed step.

Two axes, because they answer different questions:

  ADC   calibration-free: where the discriminator actually sits.  Comparing the
        thresholds on this axis says what the DAC setting does, with no MIP
        table involved.
  MIP   each set divided by its OWN MPV table.  This is the axis the
        reconstruction works on, and the one the digitiser needs, so a biased
        MPV shows up here and not in ADC.

The fit is the shape a hard discriminator on a noisy amplitude produces,

    P(x) = plateau * 0.5 * (1 + erf((x - mu) / (sqrt(2) * sigma)))

whose three parameters map onto the digitiser one for one: ``mu`` is
``Threshold``, ``sigma`` is the spread of the fast channel (``FastNoiseMIP``),
and ``plateau`` is ``TriggerEfficiency`` — the part that is not a threshold at
all but an inefficiency the chip has above it.

Usage
-----
    python -m analysis.trigger_turnon \\
        --run th230=/eos/.../TB2026CERN_run_000013/chunks \\
        --run th220=/eos/.../TB2026CERN_run_000072/chunks \\
        --run th210=/eos/.../TB2026CERN_eudaq_run_000253/chunks \\
        --events 40 --chunks 2 --outdir plots/ --write-config mappings/trigger_turnon.yml
"""

from __future__ import annotations

import argparse
import glob
import os
import sys

import numpy as np

from analysis.compare_shapers import load_tables, data_turnon, sim_turnon

# Fixed hue per threshold set, assigned in order and never cycled.
C_TH = {"th210": "#2a78d6", "th215": "#4a3aa7", "th220": "#1baf7a",
        "th230": "#eb6834"}
C_FALLBACK = ["#eda100", "#e87ba4", "#008300"]
C_GRID = "#d8d8d4"
C_INK = "#0b0b0b"
C_INK_SOFT = "#52514e"


def _run_name(chunk_path: str) -> str:
    """The run directory a chunk belongs to (<run>/chunks/chunk_NNNN.root)."""
    parent = os.path.dirname(os.path.abspath(chunk_path))
    if os.path.basename(parent) == "chunks":
        parent = os.path.dirname(parent)
    return os.path.basename(parent)


def _sim_response(path):
    """(hits per event, MIP per event) of a simulated ecal tree, unmasked hits."""
    import uproot

    arrays = uproot.open(path)["ecal"].arrays(
        ["nhit_chan", "hit_energy", "hit_ismasked"], library="np")
    n = len(arrays["nhit_chan"])
    hits = energy = 0.0
    for i in range(n):
        keep = arrays["hit_ismasked"][i] == 0
        hits += int(keep.sum())
        energy += float(arrays["hit_energy"][i][keep].sum())
    return hits / max(n, 1), energy / max(n, 1)


def turn_on_model(x, mu, sigma, plateau):
    from scipy.special import erf
    return plateau * 0.5 * (1.0 + erf((x - mu) / (np.sqrt(2.0) * sigma)))


def fit_turn_on(centres, num, den, min_entries=30):
    """Fit P(x); returns (mu, sigma, plateau) and their errors."""
    from scipy.optimize import curve_fit

    ok = den >= min_entries
    if ok.sum() < 6:
        return (np.nan,) * 3, (np.nan,) * 3
    x = centres[ok]
    eff = num[ok] / den[ok]
    # Binomial errors, floored so empty-ish bins do not dominate the fit.
    err = np.sqrt(np.maximum(eff * (1 - eff), 1e-4) / den[ok])

    plateau0 = float(np.mean(eff[x > np.percentile(x, 70)]))
    mu0 = float(x[np.argmin(np.abs(eff - 0.5 * plateau0))])
    try:
        popt, pcov = curve_fit(turn_on_model, x, eff, p0=[mu0, 0.15 * mu0, plateau0],
                               sigma=err, absolute_sigma=True, maxfev=20000,
                               bounds=([0.0, 1e-4, 0.0], [np.inf, np.inf, 1.0]))
        return tuple(popt), tuple(np.sqrt(np.diag(pcov)))
    except Exception as exc:                      # noqa: BLE001
        print(f"  [fit] failed ({exc}); falling back to the raw crossings")
        return (mu0, 0.15 * mu0, plateau0), (np.nan,) * 3


def _style(ax, xlabel, ylabel, title=None):
    ax.set_xlabel(xlabel, color=C_INK)
    ax.set_ylabel(ylabel, color=C_INK)
    if title:
        ax.set_title(title, color=C_INK, fontsize=11, pad=10)
    ax.grid(True, color=C_GRID, linewidth=0.6)
    ax.set_axisbelow(True)
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)
    for spine in ("left", "bottom"):
        ax.spines[spine].set_color(C_GRID)
    ax.tick_params(colors=C_INK_SOFT)


def per_slab_fits(centres, num_den, all_fit, override_sigma):
    """Fit the ADC-axis turn-on of every slab; flag the outliers.

    Returns ``{"fits": {slab: (mu, sigma, plateau)}, "outliers": [slab, ...]}``.
    A slab is an outlier when its discriminator level sits further from the
    all-slab value than ``override_sigma`` times the scatter of the other slabs
    -- slab 12 of the 2026 stack, at DAC 243 against 215-230, is the case this
    exists for.
    """
    num, den = num_den
    fits = {}
    for slab in range(num.shape[0]):
        if den[slab].sum() < 2000:
            continue
        try:
            par, _ = fit_turn_on(centres, num[slab], den[slab])
        except Exception:   # noqa: BLE001 -- an unfittable slab is reported, not fatal
            continue
        fits[slab] = tuple(float(v) for v in par)
    if not fits:
        return {"fits": {}, "outliers": []}
    mus = np.array([f[0] for f in fits.values()])
    med = float(np.median(mus))
    mad = float(np.median(np.abs(mus - med))) * 1.4826 or 1.0
    outliers = [s for s, f in fits.items() if abs(f[0] - med) > override_sigma * mad]
    print(f"  per slab [ADC]: {'slab':>4} {'mu':>7} {'sigma':>7} {'plateau':>8}   "
          f"(all slabs: mu = {all_fit[0]:.2f}; median of slabs {med:.2f}, scatter {mad:.2f})")
    for slab, (mu, sigma, plateau) in sorted(fits.items()):
        flag = "  <-- own discriminator" if slab in outliers else ""
        print(f"                  {slab:>4} {mu:>7.2f} {sigma:>7.2f} {plateau:>8.3f}{flag}")
    return {"fits": fits, "outliers": outliers}


def draw_per_slab(th, r, centres, outdir, tag):
    import matplotlib.pyplot as plt

    fits = r["per_slab"]["fits"]
    outliers = sorted(r["per_slab"]["outliers"])
    fig, ax = plt.subplots(figsize=(8.4, 5.4))
    fine = np.linspace(centres[0], centres[-1], 400)
    # Every highlighted curve gets its own colour AND its own dash pattern, so
    # two outlier slabs stay apart in greyscale too; the all-slab fit takes the
    # last style of the cycle.
    styles = [("#17795a", "-"), ("#2a78d6", "--"), ("#d2551f", ":")]
    for slab, (mu, sigma, plateau) in sorted(fits.items()):
        if slab in outliers:
            continue
        ax.plot(fine, turn_on_model(fine, mu, sigma, plateau), color=C_INK_SOFT,
                linewidth=1.0, alpha=0.45)
    for i, slab in enumerate(outliers):
        mu, sigma, plateau = fits[slab]
        colour, dash = styles[i % len(styles)]
        ax.plot(fine, turn_on_model(fine, mu, sigma, plateau), color=colour,
                linewidth=2.2, linestyle=dash, label=f"slab {slab}: {mu:.1f} ADC")
    mu, sigma, plateau = r["fit_adc"]
    colour, dash = styles[len(outliers) % len(styles)]
    ax.plot(fine, turn_on_model(fine, mu, sigma, plateau), color=colour,
            linewidth=2.8, linestyle=dash, label=f"all slabs: {mu:.1f} ADC")
    ax.plot([], [], color=C_INK_SOFT, linewidth=1.0, alpha=0.45, label="other slabs")
    ax.axhline(1.0, color=C_INK_SOFT, linewidth=1.0, linestyle=":")
    ax.set_xlim(0, 80)
    ax.set_ylim(-0.03, 1.08)
    _style(ax, "slow-shaper amplitude  [ADC, pedestal-subtracted]",
           "P(cell enters the event)", f"Trigger turn-on per slab — {th}, {r['source']}")
    ax.legend(frameon=False, labelcolor=C_INK, loc="lower right")
    fig.tight_layout()
    path = os.path.join(outdir, f"trigger_turnon_{th}_per_slab_{tag}.png")
    fig.savefig(path, dpi=140)
    plt.close(fig)
    return path


def main(argv=None) -> int:
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--run", action="append", default=[], required=True,
                   metavar="TH=CHUNKDIR",
                   help="Threshold set and the decoded chunk directory (or a "
                        "single chunk file), repeatable")
    p.add_argument("--events", type=int, default=40,
                   help="Acquisitions read per chunk (default: 40)")
    p.add_argument("--chunks", type=int, default=2,
                   help="Chunks read per run (default: 2)")
    p.add_argument("--calib-dir",
                   default="../siwecal-tb2026/calibration/MuonCalib_gaudi")
    p.add_argument("--outdir", default="plots")
    p.add_argument("--tag", default="allth")
    p.add_argument("--sim", action="append", default=[], metavar="TH=REAL[,SIMPLE]",
                   help="Simulated ecal trees for that threshold set. With both "
                        "chains (real and simple, same sample) the digitiser's own "
                        "turn-on is overlaid; with only the real one just the "
                        "response is taken, since the turn-on needs the simple "
                        "chain as the denominator")
    p.add_argument("--sim-label", default="simulation",
                   help="Legend entry for the simulated points")
    p.add_argument("--data-point", action="append", default=[],
                   metavar="TH=HITS,MIP,LABEL",
                   help="Measured hits/event and MIP/event of a data run, drawn "
                        "on the response plot; repeatable, one per threshold set "
                        "(e.g. th230=453,2241,run7)")
    p.add_argument("--per-slab", action="store_true",
                   help="Also fit the ADC-axis turn-on slab by slab and draw it "
                        "(trigger_turnon_<th>_per_slab_<tag>.png): the way to see "
                        "a slab that ran at a discriminator of its own")
    p.add_argument("--override-sigma", type=float, default=3.0,
                   help="--per-slab: slabs whose discriminator level is further than "
                        "this many sigma (of the per-slab scatter) from the "
                        "all-slab fit are written as slab_overrides (default: 3)")
    p.add_argument("--write-config", default="",
                   help="Write the fitted parameters as YAML for job3")
    args = p.parse_args(argv)

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    os.makedirs(args.outdir, exist_ok=True)

    bins_mip = np.arange(0.0, 5.001, 0.1)
    bins_adc = np.arange(0.0, 160.1, 2.0)
    c_mip = 0.5 * (bins_mip[:-1] + bins_mip[1:])
    c_adc = 0.5 * (bins_adc[:-1] + bins_adc[1:])

    results = {}
    for i, spec in enumerate(args.run):
        if "=" not in spec:
            sys.exit(f"ERROR: --run expects TH=CHUNKDIR, got '{spec}'")
        th, source = spec.split("=", 1)
        chunks = (sorted(glob.glob(os.path.join(source, "chunk_*.root")))[:args.chunks]
                  if os.path.isdir(source) else [source])
        if not chunks:
            sys.exit(f"ERROR: no chunks under {source}")

        print(f"\n=== {th}: {len(chunks)} chunk(s) from {source}")
        mpv, ped = load_tables(args.calib_dir, th)
        out = data_turnon(chunks, args.events, mpv, ped, bins_mip, bins_adc)
        results[th] = {
            "mip": out["mip"], "adc": out["adc"], "n": out["n"],
            "mpv": float(np.nanmedian(mpv)),
            "color": C_TH.get(th, C_FALLBACK[i % len(C_FALLBACK)]),
            # chunks live in <run>/chunks/, so the run name is two levels up.
            "source": _run_name(chunks[0]),
        }
        for axis, centres in (("mip", c_mip), ("adc", c_adc)):
            num, den = results[th][axis]
            par, err = fit_turn_on(centres, num, den)
            results[th][f"fit_{axis}"] = par
            results[th][f"err_{axis}"] = err
            unit = "MIP" if axis == "mip" else "ADC"
            print(f"  turn-on [{unit}]: mu = {par[0]:.2f}, sigma = {par[1]:.2f}, "
                  f"plateau = {par[2]:.3f}")
        if args.per_slab and out["adc_per_slab"] is not None:
            results[th]["per_slab"] = per_slab_fits(c_adc, out["adc_per_slab"],
                                                     results[th]["fit_adc"], args.override_sigma)

    # ------------------------------------------------- simulated counterpart
    sims = {}
    for spec in args.sim:
        if "=" not in spec:
            sys.exit(f"ERROR: --sim expects TH=REAL[,SIMPLE], got '{spec}'")
        th, paths = spec.split("=", 1)
        real_path, _, simple_path = paths.partition(",")
        print(f"\n=== {th}: simulation")
        sims[th] = {"mip": None, "fit": None,
                    "response": _sim_response(real_path)}
        # The turn-on is P(kept by the real chain | simple-chain amplitude), so it
        # needs the simple chain as the denominator.  The same file twice would
        # give an efficiency of 1 everywhere, which is not a measurement.
        if simple_path and simple_path != real_path:
            num, den = sim_turnon(real_path, simple_path, bins_mip)
            par, _ = fit_turn_on(c_mip, num, den)
            sims[th].update({"mip": (num, den), "fit": par})
            print(f"  turn-on [MIP]: mu = {par[0]:.2f}, sigma = {par[1]:.2f}, "
                  f"plateau = {par[2]:.3f}")
        else:
            print("  turn-on: skipped (no simple chain given)")
        print(f"  response     : {sims[th]['response'][0]:.1f} hits/ev, "
              f"{sims[th]['response'][1]:.0f} MIP/ev")

    # ------------------------------------------------------------------ plots
    written = []
    for axis, centres, xlabel, xmax in (
            ("adc", c_adc, "slow-shaper amplitude  [ADC, pedestal-subtracted]", 120),
            ("mip", c_mip, "slow-shaper amplitude  [MIP, that set's own table]", 3.5)):
        fig, ax = plt.subplots(figsize=(8.4, 5.4))
        for th, r in results.items():
            num, den = r[axis]
            eff = np.where(den >= 30, num / np.maximum(den, 1), np.nan)
            err = np.where(den >= 30,
                           np.sqrt(np.maximum(eff * (1 - eff), 0) / np.maximum(den, 1)),
                           np.nan)
            mu, sigma, plateau = r[f"fit_{axis}"]
            ax.errorbar(centres, eff, yerr=err, color=r["color"], linewidth=0,
                        elinewidth=1.2, marker="o", markersize=3.5,
                        label=f"{th} — {r['source']}")
            fine = np.linspace(centres[0], centres[-1], 400)
            ax.plot(fine, turn_on_model(fine, mu, sigma, plateau), color=r["color"],
                    linewidth=1.8, alpha=0.85)
        if axis == "mip":
            for th, sim in sims.items():
                if sim["mip"] is None:
                    continue
                num, den = sim["mip"]
                eff = np.where(den >= 30, num / np.maximum(den, 1), np.nan)
                ax.plot(centres, eff, color=results[th]["color"], linewidth=1.6,
                        linestyle="--", label=f"{th} — simulation")
        ax.axhline(1.0, color=C_INK_SOFT, linewidth=1.0, linestyle=":")
        ax.set_xlim(0, xmax)
        ax.set_ylim(-0.03, 1.08)
        unit = "ADC" if axis == "adc" else "MIP"
        _style(ax, xlabel, "P(cell enters the event)",
               f"Trigger turn-on per threshold set — {unit} axis")
        ax.legend(frameon=False, labelcolor=C_INK, loc="lower right")
        fig.tight_layout()
        path = os.path.join(args.outdir, f"trigger_turnon_{axis}_{args.tag}.png")
        fig.savefig(path, dpi=140)
        plt.close(fig)
        written.append(path)

    if sims:
        # Ordered by the DISCRIMINATOR level in ADC -- the calibration-free
        # quantity -- but plotted per set, not against it: the configured MIP
        # threshold is each set's ADC level divided by its own MPV table, so the
        # two orderings differ and an ADC x-axis would suggest a monotonicity the
        # numbers do not have.
        order = sorted(sims, key=lambda t: results[t]["fit_adc"][0])
        x = np.arange(len(order))
        fig, (ax, axh) = plt.subplots(2, 1, figsize=(8.4, 6.8), sharex=True)
        ax.plot(x, [sims[th]["response"][1] for th in order], color=C_TH["th230"],
                linewidth=2.0, marker="o", markersize=7, label=args.sim_label)
        axh.plot(x, [sims[th]["response"][0] for th in order], color=C_TH["th230"],
                 linewidth=2.0, marker="o", markersize=7)
        for i, th in enumerate(order):
            mu, sigma, plateau = results[th]["fit_mip"]
            if not np.isfinite(mu) or not np.isfinite(plateau):
                continue
            ax.annotate(f"threshold {mu:.2f} MIP\nefficiency {plateau:.2f}",
                        (i, sims[th]["response"][1]), textcoords="offset points",
                        xytext=(0, 12), ha="center", fontsize=8.5, color=C_INK_SOFT)
        labelled = False
        for spec in args.data_point:
            dth, values = spec.split("=", 1)
            parts = values.split(",")
            hits, mip = float(parts[0]), float(parts[1])
            note = parts[2] if len(parts) > 2 else ""
            if dth not in order:
                continue
            i = order.index(dth)
            ax.plot([i], [mip], color=C_TH["th210"], marker="s", markersize=9,
                    linestyle="none",
                    label=None if labelled else "test beam")
            axh.plot([i], [hits], color=C_TH["th210"], marker="s", markersize=9,
                     linestyle="none")
            if note:
                ax.annotate(note, (i, mip), textcoords="offset points",
                            xytext=(14, -4), ha="left", fontsize=8.5,
                            color=C_TH["th210"])
            labelled = True
        labels = [f"{th}\n{results[th]['fit_adc'][0]:.1f} ADC"
                  if np.isfinite(results[th]["fit_adc"][0]) else th
                  for th in order]
        axh.set_xticks(x)
        axh.set_xticklabels(labels)
        ax.set_ylim(bottom=0)
        axh.set_ylim(bottom=0)
        _style(ax, "", "event energy  [MIP]",
               "The same sample digitised with each threshold set")
        _style(axh, "threshold set (measured discriminator level)",
               "hits per event")
        ax.legend(frameon=False, labelcolor=C_INK, loc="lower left")
        fig.tight_layout()
        path = os.path.join(args.outdir, f"trigger_response_{args.tag}.png")
        fig.savefig(path, dpi=140)
        plt.close(fig)
        written.append(path)

    # The DAC law: if the trigger DAC only moves the discriminator and nothing
    # else, these points lie on a line.  Three points is not a lot, but it is the
    # independent check on the same statement the pedestals and the LG->HG anchor
    # make -- and it turns the digitiser's threshold input into a measured law
    # rather than one number per calibration set.
    dacs = []
    for th in results:
        digits = "".join(ch for ch in th if ch.isdigit())
        if digits and np.isfinite(results[th]["fit_adc"][0]):
            dacs.append((float(digits), results[th]["fit_adc"][0],
                         results[th]["err_adc"][0], th))
    if len(dacs) >= 2:
        dacs.sort()
        x = np.array([d[0] for d in dacs])
        y = np.array([d[1] for d in dacs])
        err = np.array([d[2] if np.isfinite(d[2]) else 0.0 for d in dacs])
        slope, intercept = np.polyfit(x, y, 1)
        fig, ax = plt.subplots(figsize=(8.0, 5.0))
        for xi, yi, ei, th in dacs:
            ax.errorbar([xi], [yi], yerr=[ei], color=results[th]["color"],
                        marker="o", markersize=8, linestyle="none", capsize=0,
                        label=th)
        fine = np.linspace(x.min() - 5, x.max() + 5, 100)
        ax.plot(fine, slope * fine + intercept, color=C_INK_SOFT, linewidth=1.5,
                linestyle="--",
                label=f"{slope:.3f} ADC per DAC unit")
        _style(ax, "trigger threshold DAC", "discriminator level  [ADC]",
               "What the DAC setting is worth in ADC")
        ax.legend(frameon=False, labelcolor=C_INK)
        fig.tight_layout()
        path = os.path.join(args.outdir, f"threshold_vs_dac_{args.tag}.png")
        fig.savefig(path, dpi=140)
        plt.close(fig)
        written.append(path)
        residual = y - (slope * x + intercept)
        print(f"\n[dac ] {slope:.3f} ADC per DAC unit, intercept {intercept:.1f} ADC; "
              f"residuals " + ", ".join(f"{th}: {r:+.2f}" for (_, _, _, th), r
                                        in zip(dacs, residual)))

    # ---------------------------------------------------------------- summary
    lines = [f"# Trigger turn-on per threshold set "
             f"({args.chunks} chunk(s) x {args.events} acquisitions each)",
             "",
             f"{'th':<8}{'run':<26}{'MPV':>8}{'mu[ADC]':>10}{'sig[ADC]':>10}"
             f"{'mu[MIP]':>10}{'sig[MIP]':>10}{'plateau':>10}{'mu/MPV':>9}"]
    for th, r in results.items():
        a, m = r["fit_adc"], r["fit_mip"]
        lines.append(f"{th:<8}{r['source']:<26}{r['mpv']:>8.2f}{a[0]:>10.2f}"
                     f"{a[1]:>10.2f}{m[0]:>10.2f}{m[1]:>10.2f}{m[2]:>10.3f}"
                     f"{a[0] / r['mpv']:>9.2f}")
    summary = "\n".join(lines)
    print("\n" + summary)
    summary_path = os.path.join(args.outdir, f"trigger_turnon_{args.tag}.txt")
    with open(summary_path, "w") as fh:
        fh.write(summary + "\n")
    written.append(summary_path)

    if args.per_slab:
        for th, r in results.items():
            if r.get("per_slab"):
                written.append(draw_per_slab(th, r, c_adc, args.outdir, args.tag))

    if args.write_config:
        import yaml
        config = {
            "_comment": "Measured by analysis/trigger_turnon.py from hitbit_high "
                        "in the decoded chunks. mu/sigma are in MIP on that "
                        "threshold set's own MPV table; RealDigitizer reads them "
                        "as Threshold / FastNoiseMIP / TriggerEfficiency.",
            "thresholds": {
                th: {"threshold_mip": round(float(r["fit_mip"][0]), 4),
                     "sigma_mip": round(float(r["fit_mip"][1]), 4),
                     "efficiency": round(float(r["fit_mip"][2]), 4),
                     "threshold_adc": round(float(r["fit_adc"][0]), 2),
                     "source_run": r["source"]}
                for th, r in results.items()},
        }
        # Slabs at a discriminator of their own, in ADC (their gain is normal, so
        # job3 turns them into MIP with the SET's adc_per_mip).
        for th, r in results.items():
            ps = r.get("per_slab")
            if not ps or not ps["outliers"]:
                continue
            config["thresholds"][th]["slab_overrides"] = {
                int(slab): {"threshold_adc": round(ps["fits"][slab][0], 2),
                            "sigma_adc": round(ps["fits"][slab][1], 2),
                            "efficiency": round(ps["fits"][slab][2], 4),
                            "turnon_source_run": r["source"]}
                for slab in ps["outliers"]}
        os.makedirs(os.path.dirname(os.path.abspath(args.write_config)), exist_ok=True)
        with open(args.write_config, "w") as fh:
            yaml.safe_dump(config, fh, sort_keys=False)
        written.append(args.write_config)

    print("\nWritten:")
    for path in written:
        print(f"  {path}")
    return 0


if __name__ == "__main__":
    sys.exit(main())

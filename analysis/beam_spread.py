#!/usr/bin/env python
"""
Beam spot size per group and per energy, measured from the shower barycentre (bar_x, bar_y) of data and of the
final_v5 simulation, for the beam sigmas of the next ddsim campaign (simulation/run_script/launch_beam_final.sh).

    python analysis/beam_spread.py [--out analysis/beam_spread.txt]
    # closure on the new campaign: the beam it returns must be the one it was simulated with
    python analysis/beam_spread.py --groups analysis/final_v6_groups.txt --sim .../final_v6 \
        --sim-sigma analysis/beam_spread.txt --out final_v6_beam_spread.txt

The barycentre spread is the beam spot convolved with the shower-barycentre resolution (and the acceptance, which
data and simulation share because the sim sample sits at the group's position). The simulated beam is Gaussian with
known sigma (13.75 x 8.25 mm in final_v5, per energy from --sim-sigma otherwise), so the data beam follows by quadrature:
    sigma_beam^2 = sigma_sim_beam^2 + s_data^2 - s_sim^2      (s = robust width of bar_x or bar_y, IQR / 1.349)
Selection as in final_panels.py: shower-like events of the valtrees (nhit_chan > 0.5 x p90).
Per energy: the median over the groups at that energy (positions near the edge truncate the distribution; the median
keeps them from pulling the value).
"""
import argparse
import glob
import os
import sys

import numpy as np
import uproot

_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(_HERE))
from analysis.final_panels import GROUPS_FILE, RECO, SIM, load_groups, showerlike  # noqa: E402

SIM_SIGMA = (13.75, 8.25)  # beam sigma (x, y) [mm] of the final_v5 ddsim samples


def read_sigma_table(path):
    """Per-energy beam sigma {E: (sx, sy)} from the 'per energy' block of a beam_spread.txt."""
    out, on = {}, False
    for line in open(path):
        if line.startswith("# per energy"):
            on = True
        elif on and not line.startswith("#") and line.split():
            E, _, sx, sy = line.split()[:4]
            out[float(E)] = (float(sx), float(sy))
    return out


def sim_valtree(sim_dir, g, kind):
    v = glob.glob(f"{sim_dir}/pid/{g['name']}_{kind}/*.valtree.root")
    return v[0] if v else None


def barycentres(paths):
    cols = {b: [] for b in ("nhit_chan", "bar_x", "bar_y")}
    for p in paths:
        for a in uproot.iterate(f"{p}:ecal", list(cols), step_size="400 MB", library="np"):
            for b in cols:
                cols[b].append(a[b])
    c = {b: np.concatenate(v).astype(float) for b, v in cols.items()}
    k = showerlike(c["nhit_chan"])
    return c["bar_x"][k], c["bar_y"][k]


def width(v):
    q25, q75 = np.percentile(v, [25, 75])
    return (q75 - q25) / 1.349


def beam_sigma(s_data, s_sim, s_beam_sim):
    v = s_beam_sim ** 2 + s_data ** 2 - s_sim ** 2
    return np.sqrt(v) if v > 0 else 0.0


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", default=os.path.join(_HERE, "beam_spread.txt"))
    ap.add_argument("--groups", default=GROUPS_FILE)
    ap.add_argument("--sim", default=SIM, help="simulation campaign dir (pid/<group>_digi/*.valtree.root)")
    ap.add_argument("--sim-sigma", default=None,
                    help="beam_spread.txt the simulation was generated with (default: 13.75 x 8.25 for every energy)")
    a = ap.parse_args(argv)
    sim_sigma = read_sigma_table(a.sim_sigma) if a.sim_sigma else {}

    rows, per_E = [], {}
    for g in load_groups(a.groups):
        runs = [r for r in g["runs"] if os.path.exists(f"{RECO}/{r}/ecal_{r}.valtree.root")]
        sim = sim_valtree(a.sim, g, "digi")
        if not runs or sim is None:
            print(f"[{g['name']}] skipped (data {len(runs)} runs, sim {sim})")
            continue
        dx, dy = barycentres([f"{RECO}/{r}/ecal_{r}.valtree.root" for r in runs])
        sx, sy = barycentres([sim])
        w = (width(dx), width(dy), width(sx), width(sy))
        ss = sim_sigma.get(float(g["E"]), SIM_SIGMA)
        bx, by = beam_sigma(w[0], w[2], ss[0]), beam_sigma(w[1], w[3], ss[1])
        rows.append((g["name"], len(dx), len(sx), *w, bx, by))
        per_E.setdefault(float(g["E"]), []).append((bx, by))
        print(f"{g['name']:22s} data {w[0]:5.1f} {w[1]:5.1f}  sim {w[2]:5.1f} {w[3]:5.1f}  -> beam {bx:5.1f} {by:5.1f}")

    with open(a.out, "w") as f:
        f.write(f"# analysis/beam_spread.py: barycentre width (IQR/1.349) [mm], shower-like events; sim {a.sim}, "
                f"sim beam sigma {a.sim_sigma or SIM_SIGMA}; beam = sqrt(sim_beam^2 + data^2 - sim^2)\n")
        f.write(f"# {'group':22s} {'Ndata':>7s} {'Nsim':>5s} {'dat_x':>6s} {'dat_y':>6s} {'sim_x':>6s} {'sim_y':>6s} "
                f"{'beam_x':>6s} {'beam_y':>6s}\n")
        for r in rows:
            f.write(f"  {r[0]:22s} {r[1]:7d} {r[2]:5d} " + " ".join(f"{v:6.2f}" for v in r[3:]) + "\n")
        f.write("#\n# per energy (median over groups) -> launch_beam_final.sh\n# E[GeV] ngroups sigma_x sigma_y\n")
        for E in sorted(per_E):
            v = np.array(per_E[E])
            f.write(f"  {E:5g} {len(v):3d} {np.median(v[:, 0]):6.2f} {np.median(v[:, 1]):6.2f}\n")
    print(open(a.out).read().split("# per energy")[1])


if __name__ == "__main__":
    main()

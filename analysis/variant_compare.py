#!/usr/bin/env python
"""
Simulation variants against data, from valtrees: shower-like events (nhit > 0.5 x p90), Gaussian-core nhit and energy
[MIP], median energy per hit, mean hits in the core (r < 10 mm), in the halo (15-60 mm, layers 3-8) and in layers 0-2,
median z barycentre and Moliere radius; every value as sim / data.

    python analysis/variant_compare.py <out.txt> <data valtree> <label>=<sim valtree> [<label>=<sim valtree> ...]
"""
import sys

import numpy as np
import uproot

sys.path.insert(0, __file__.rsplit("/analysis/", 1)[0])
sys.path.insert(0, __file__.rsplit("/siwecal_k4sim/", 1)[0] + "/siwecal-tb2026")
from analysis.compare_event_level import gauss_core  # noqa: E402


def stats(path, nmax=30000):
    a = uproot.open(path)["ecal"].arrays(["nhit_chan", "sum_energy", "e_over_nhit", "zbary", "moliere", "bar_x", "bar_y",
                                          "hit_x", "hit_y", "hit_slab"], entry_stop=nmax, library="np")
    k = a["nhit_chan"] > 0.5 * np.percentile(a["nhit_chan"], 90)
    core = halo = front = 0.0
    idx = np.nonzero(k)[0][:6000]
    for i in idx:
        r = np.hypot(a["hit_x"][i] - a["bar_x"][i], a["hit_y"][i] - a["bar_y"][i]); s = a["hit_slab"][i]
        core += (r < 10).sum(); halo += ((r > 15) & (r < 60) & (s >= 3) & (s <= 8)).sum(); front += (s <= 2).sum()
    n = len(idx)
    return {"nhit": gauss_core(a["nhit_chan"][k].astype(float))[0], "E": gauss_core(a["sum_energy"][k].astype(float))[0],
            "E/hit": np.median(a["e_over_nhit"][k]), "core": core / n, "halo": halo / n, "L0-2": front / n,
            "zbary": np.median(a["zbary"][k]), "moliere": np.median(a["moliere"][k])}


def main():
    out, data = sys.argv[1], sys.argv[2]
    d = stats(data)
    keys = list(d)
    lines = [f"data: {data}", f"{'':<14}" + "".join(f"{k:>9}" for k in keys), f"{'data':<14}" + "".join(f"{d[k]:9.2f}" for k in keys)]
    for arg in sys.argv[3:]:
        lab, p = arg.split("=", 1)
        s = stats(p)
        lines.append(f"{lab:<14}" + "".join(f"{s[k] / d[k]:9.3f}" for k in keys))
        print(lines[-1], flush=True)
    open(out, "w").write("\n".join(lines) + "\n")
    print("\n".join(lines[:3]))


if __name__ == "__main__":
    main()

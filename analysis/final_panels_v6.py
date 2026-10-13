#!/usr/bin/env python
"""
The full 15-variable shower panel (final_panels.py's) for final_v6 groups, before the whole th210 scan exists in
final_v6: the MIP -> GeV line is final_v5's (final_panels_select.CAL), the rest is final_panels.py unchanged.

    FINAL_SIM=.../final_v6 python analysis/final_panels_v6.py <outdir> --groups th210_S1_e20,th210_S1_e74

Group names from analysis/final_v6_groups.txt (FINAL_GROUPS overrides). Writes <outdir>/panel_<group>.png and
<outdir>/panels_v6_summary.txt (Gaussian-core mu and sigma/mu, medians, sim / data).
"""
import argparse
import os
import sys

import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(_HERE))
os.environ.setdefault("FINAL_GROUPS", os.path.join(_HERE, "final_v6_groups.txt"))
from analysis.final_panels import KINDS, PANELS, group_samples, load_groups, panel  # noqa: E402
from analysis.final_panels_select import CAL  # noqa: E402


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("outdir")
    ap.add_argument("--groups", required=True, help="Comma-separated final_v6 group names")
    ap.add_argument("--no-cache", action="store_true")
    a = ap.parse_args(argv)
    os.makedirs(a.outdir, exist_ok=True)
    want = a.groups.split(",")
    groups = [g for g in load_groups(os.environ["FINAL_GROUPS"]) if g["name"] in want]
    missing = set(want) - {g["name"] for g in groups}
    if missing:
        sys.exit(f"unknown groups: {sorted(missing)}")
    lines = ["MIP -> GeV: final_v5's line per sample kind (final_panels_select.CAL)"]
    for g in groups:
        s = group_samples(g, a.outdir, not a.no_cache)
        for k in KINDS:
            s[k]["gev"] = (s[k]["sum_w"] - CAL[k][1]) / CAL[k][0]
        st = panel(g, s, f"{a.outdir}/panel_{g['name']}.png")
        lines.append(f"\n== {g['name']}  N data/nodigi/digi = "
                     f"{st[('nhit_chan', 'data')][2]}/{st[('nhit_chan', 'nodigi')][2]}/{st[('nhit_chan', 'digi')][2]}")
        for v, _ in PANELS:
            m = {k: st.get((v, k), (np.nan, np.nan, 0)) for k in KINDS}
            r = (m["digi"][0] / m["data"][0]) if m["data"][0] else np.nan
            lines.append(f"   {v:<18}" + "".join(f"{m[k][0]:10.4g}" for k in KINDS) + f"   digi/data {r:.3f}")
    open(f"{a.outdir}/panels_v6_summary.txt", "w").write("\n".join(lines) + "\n")
    print("\n".join(lines))


if __name__ == "__main__":
    main()

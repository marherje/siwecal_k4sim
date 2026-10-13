"""One representative 20 GeV e+ shower from data (eudaq 291, th210) and one from the digitised simulation (final_v5),
drawn with the event viewer's own 3D scene plus the Molière cylinder and the shower axis (barycentre).

Representative = shower-like (nhit > 0.5 x p90) and closest to the sample's medians of energy, Molière radius and
longitudinal and transverse barycentre (distance in units of each variable's interquartile range), so both
showers sit at the centre of their beam spot.
"""
import os
import sys

import numpy as np
import uproot
import plotly.graph_objects as go

sys.path.insert(0, "/afs/cern.ch/user/m/marquezh/public/siwecal-tb2026")
from event_viewer.config import ViewerConfig  # noqa: E402
from event_viewer.controller import ViewerController  # noqa: E402

B = "/eos/experiment/drdcalo/siw-ecal/TB2026-06"
SAMPLES = [
    ("data", "Data, eudaq 291", f"{B}/Reconstructed_final/TB2026CERN_eudaq_run_000291/ecal_TB2026CERN_eudaq_run_000291.valtree.root"),
    ("sim", "Simulation (digitised, final_v5)", f"{B}/Simulation/Processed/adc_vs_tb/final_v5/pid/th210_e20_26_-52_digi/ecal_th210_e20_26_-52.valtree.root"),
]
OUT = os.path.dirname(os.path.abspath(__file__))
VARS = ["energy", "moliere", "zbary", "bar_x", "bar_y"]
HALF = 90.0


def representative(path):
    a = uproot.open(path)["ecal"].arrays(["nhit"] + VARS, library="np")
    ok = (a["nhit"] > 0.5 * np.percentile(a["nhit"], 90)) & np.all([np.isfinite(a[v]) for v in VARS], axis=0)
    d = np.zeros(len(a["nhit"]))
    for v in VARS:
        x = a[v]
        q1, med, q3 = np.percentile(x[ok], [25, 50, 75])
        d += ((x - med) / max(q3 - q1, 1e-9)) ** 2
    d[~ok] = np.inf
    i = int(np.argmin(d))
    return i, {v: float(np.median(a[v][ok])) for v in VARS}, {v: float(a[v][i]) for v in VARS + ["nhit"]}


ctl = ViewerController(ViewerConfig())
for tag, label, path in SAMPLES:
    idx, med, val = representative(path)
    ds = ctl.dataset(path)
    ev = ds.get_event(idx)
    fig = ctl.scene3d.event_figure(ev, show_moliere=True, show_axis=True)
    # The barycentre itself, as a clear marker: (bar_x, bar_y) and z_bary mapped from layers to mm
    m = ev.metrics
    slab_z = ctl.detector.slab_z_mm
    zb = float(np.interp(m["zbary"], np.arange(slab_z.size), slab_z))
    fig.add_trace(go.Scatter3d(x=[m["bar_x"]], y=[m["bar_y"]], z=[zb], mode="markers",
                               marker=dict(size=9, color="#d62728", symbol="diamond", line=dict(color="white", width=2)),
                               name=f"barycentre ({m['bar_x']:.0f}, {m['bar_y']:.0f}) mm, z_bary {m['zbary']:.2f} layers"))
    # Same framing on both: +-90 mm around the shower axis in x and y, the full depth in z
    fig.update_layout(scene=dict(xaxis=dict(range=[m["bar_x"] - HALF, m["bar_x"] + HALF]),
                                 yaxis=dict(range=[m["bar_y"] - HALF, m["bar_y"] + HALF]),
                                 aspectmode="data"))
    fig.update_layout(
        showlegend=True, legend=dict(x=0.01, y=0.99, bgcolor="rgba(255,255,255,0.8)"),
        title=dict(text=f"20 GeV e⁺, th210 — {label}, event {idx}: {int(val['nhit'])} hits, "
                        f"E = {val['energy']:.0f} MIP, Molière = {val['moliere']:.1f} mm, "
                        f"z<sub>bary</sub> = {val['zbary']:.2f} layers", x=0.5, y=0.97, font=dict(size=15)),
        margin=dict(l=0, r=0, t=50, b=0))
    fig.write_html(f"{OUT}/shower20_{tag}_ev{idx}.html", include_plotlyjs="cdn")
    fig.write_image(f"{OUT}/shower20_{tag}_ev{idx}.png", width=1400, height=1000, scale=1)
    print(f"{tag}: event {idx}  " + "  ".join(f"{v} {val[v]:.2f} (median {med[v]:.2f})" for v in VARS))

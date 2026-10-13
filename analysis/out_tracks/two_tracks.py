"""Find run-60 events with two ACTS tracks and render the event viewer's 3D scene to PNG + HTML."""
import os
import sys

import numpy as np
import uproot
import awkward as ak

sys.path.insert(0, "/afs/cern.ch/user/m/marquezh/public/siwecal-tb2026")
from event_viewer.config import ViewerConfig  # noqa: E402
from event_viewer.controller import ViewerController  # noqa: E402

PATH = "/eos/experiment/drdcalo/siw-ecal/TB2026-06/Reconstruction/TB2026CERN_run_000060/ecal_TB2026CERN_run_000060.edm4hep.root"
OUT = sys.argv[1]
os.makedirs(OUT, exist_ok=True)

n = ak.to_numpy(ak.num(uproot.open(PATH)["events"]["ACTSTracks/ACTSTracks.ndf"].array()))
cand = np.where(n == 2)[0]
print(f"{len(cand)} frames with 2 tracks")

ctl = ViewerController(ViewerConfig())
ds = ctl.dataset(PATH)
print("dataset events:", ds.n_events)

scored = []
for idx in cand[:200]:
    ev = ds.get_event(int(idx), want_tracks=True)
    trs = [t for t in (ev.tracks or []) if len(t.get("points") or []) >= 2]
    if len(trs) != 2:
        continue
    npts = [len(t["points"]) for t in trs]
    p0, p1 = (np.asarray(t["points"], float) for t in trs)
    m = min(len(p0), len(p1))
    sep = float(np.median(np.hypot(p0[:m, 0] - p1[:m, 0], p0[:m, 1] - p1[:m, 1])))
    scored.append((min(npts), sep, int(idx), ev.n_hits, [(t.get("chi2"), t.get("ndf")) for t in trs]))

# long tracks first, then well separated (but both inside the detector view)
good = sorted([s for s in scored if s[1] > 25], key=lambda s: (-s[0], -s[1]))
for s in good[:10]:
    print(s)
best = good[0][2] if good else scored[0][2]
for rank, idx in enumerate([s[2] for s in good[:3]] or [best]):
    ev = ds.get_event(idx, want_tracks=True)
    fig = ctl.scene3d.event_figure(ev, show_tracks=True)
    fig.update_layout(showlegend=True, legend=dict(x=0.01, y=0.99, bgcolor="rgba(255,255,255,0.8)"),
                      title=dict(text=f"Run 60 (th220 muons), frame {idx}: two ACTS tracks, {ev.n_hits} hits",
                                 x=0.5, y=0.97), margin=dict(l=0, r=0, t=40, b=0))
    fig.write_html(f"{OUT}/run60_frame{idx}.html", include_plotlyjs="cdn")
    fig.write_image(f"{OUT}/run60_frame{idx}.png", width=1400, height=1000, scale=1)
    print("wrote", idx)

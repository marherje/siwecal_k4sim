import argparse
import math
import os
import sys
from pathlib import Path
from collections import defaultdict

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import ROOT

_PLANE_BIT_OFFSET = {"SiTarget": 20, "MTCSciFi": 22, "MTCScint": 22}
_PLANE_BIT_WIDTH  = {"SiTarget": 1,  "MTCSciFi": 2,  "MTCScint": 2}

# Detectors whose hit positions are used to build the unified muon truth path.
# MTC scintillator is intentionally excluded because its CartesianGridXY pads
# are coarse (~50 mm steps) and produce a stair-step pattern that would
# dominate the otherwise smooth SciFi U/V-recovered path.  Scint hits are
# still drawn as scatter so the muon's footprint there is visible.
_TRUTH_PATH_DETS = ("SiTarget", "SiPad", "MTCSciFi")


def _hit_plane(cell_id, det_name):
    """Decode the 'plane' field from a SimCalorimeterHit cellID, or 0 if N/A."""
    off = _PLANE_BIT_OFFSET.get(det_name)
    if off is None:
        return 0  # SiPad has a single readout, no plane field
    mask = (1 << _PLANE_BIT_WIDTH[det_name]) - 1
    return int((cell_id >> off) & mask)


def _open_podio_reader(edm4hep_file):
    try:
        from podio import reading
        return reading.get_reader(edm4hep_file)
    except Exception:
        import podio
        return podio.root_io.Reader(edm4hep_file)




def read_contribs_by_pdg(edm4hep_file, det_name):
    """Read SimCalorimeterHit positions for one detector, keyed by per-hit PDG.

    Returns dict {window_id: {pdg: [(x, y, z, t, E), ...]}}

    Note on positions: the DDG4 -> EDM4HEP mapper in key4hep-2026-02-01 does not
    transfer CaloHitContribution.stepPosition, so we use hit.getPosition() (per-hit)
    instead.  For SciFi U/V strips that is (strip_center_x, energy-weighted step-y,
    layer-z); for SiTarget/SiPad it is the pad/strip center.

    Per-hit PDG: each hit aggregates many contributions; we tag the hit by the
    PDG of its highest-energy contribution.  Contribution PDGs come from the
    parallel '{det}ContribPDGs' frame parameter (vector<int>, same length and
    order as the *_Contributions collection).
    """
    reader = _open_podio_reader(edm4hep_file)
    out = defaultdict(list)
    for iframe, frame in enumerate(reader.get("events")):
        if iframe > 1: break
        try:
            hits = frame.get(f"{det_name}HitsWindowed")
            contribs = frame.get(f"{det_name}HitsWindowed_Contributions")
        except Exception:
            continue
        if hits is None or contribs is None:
            continue


        for hit in hits:

            pos = hit.getPosition()  # mm
            try:
                t = float(next(iter(hit.getContributions())).getTime())
            except StopIteration:
                t = 0.0
            E = float(hit.getEnergy())
            plane = _hit_plane(int(hit.getCellID()), det_name)
            x, y, z = float(pos.x), float(pos.y), float(pos.z)
            print(iframe, x,y,z, E, plane)
            out["iframe"].append(iframe)
            out["x"].append(x)
            out["y"].append(y)
            out["z"].append(z)
            out["E"].append(E)
            out["plane"].append(plane)
    return out
out = pd.DataFrame(read_contribs_by_pdg("timewindows.edm4hep.root", "MTCScint"))
event = out.loc[out["iframe"] == 1]

fig = plt.figure()
ax = fig.add_subplot(projection='3d')

ax.scatter(event["z"], event["x"], event["y"])
ax.set_xlabel("z")
ax.set_ylabel("x")
ax.set_zlabel("y")
plt.savefig("event_1.pdf")
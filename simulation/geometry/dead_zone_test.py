#!/usr/bin/env python3
"""Prove that the dead regions really kill hits, by shooting muons through them.

The static checks in analysis/tests/test_wafer_geometry.py show that geometry and
readout agree. They cannot show that Geant4 drops energy deposited in a gap --
that depends on which volume carries the sensitive detector, which is exactly the
thing that is easy to get wrong. This script closes that hole.

Method: fire N muons along z at a fixed (x, y) and measure the plane efficiency,
i.e. the fraction of (event, layer) pairs that produced at least one hit. A
minimum-ionising muon crossing active silicon fires every plane, so aiming at a
pad must give 100%; aiming at a dead region must make it collapse. Anything in
between means the gap is modelled in the geometry but not respected by the
readout, or the other way round.

The dead region between two sensors is their two 0.61 mm rims, so the scan below
walks in from a pad centre to the middle of the cross and the efficiency has to
collapse between 0.7 mm and 0.5 mm off the middle -- that is where the rim
starts. See the table printed by the script for the measured numbers.

Read the 0.7 mm row as a transition, not as a failure: the aim point is on pad
15, but only 0.09 mm from the rim, and a muon wanders further than that by
multiple scattering before it reaches the back of the stack. It comes out around
80%, between the 100% of the pad centre and the ~30% in the middle of the cross.

The pad-boundary row is the control that makes the test meaningful: inside a
sensor the pads butt against each other, so a muon on a pad boundary still fires
every plane, it just splits between the two neighbours. Only the rows inside the
rim lose hits. Without that control you cannot tell a correct dead region from a
readout that silently drops hits at every pad boundary.

Usage:
    source /cvmfs/sw.hsf.org/key4hep/setup.sh -r 2026-02-01
    export LD_LIBRARY_PATH=$PWD/install/lib64:$PWD/install/lib:$LD_LIBRARY_PATH
    python simulation/geometry/dead_zone_test.py                  # default scan
    python simulation/geometry/dead_zone_test.py --x 0 --y -58.8  # single point
    python simulation/geometry/dead_zone_test.py --compact simulation/geometry/other.xml
"""

from __future__ import annotations

import argparse
import collections
import os
import subprocess
import sys
import tempfile

REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
DEFAULT_COMPACT = os.path.join(REPO_ROOT, "simulation", "geometry", "SND_compact.xml")
READOUT = "SiPadHits"

# (x [mm], y [mm], label) -- the default scan, see the table in the docstring.
# Pad centres for the 1-ASU layout: pad i sits at
#     -89.09 + (i % 16 + 0.5) * 5.53 + 89.70 * (i // 16)
# i.e. the 5.53 mm pixel pitch inside a sensor, and 89.70 mm from one sensor to
# the next. The dead cross is at x = 0, 0.61 mm wide either side.
ROW_Y = -58.675  # pad row 5, kept away from any boundary in y
DEFAULT_TARGETS = [
    (-58.675, ROW_Y, "pad centre (5,5)"),
    (-3.375, ROW_Y, "pad centre (15,5), last before the dead cross"),
    (-6.140, ROW_Y, "boundary between pads 14 and 15 (inside a sensor)"),
    (-0.700, ROW_Y, "0.7 mm from the middle: still on pad 15, past the rim"),
    (-0.500, ROW_Y, "0.5 mm from the middle: inside the 0.61 mm rim"),
    (0.000, ROW_Y, "middle of the dead cross between two sensors"),
]


def simulate(compact: str, x: float, y: float, events: int, out: str) -> None:
    cmd = [
        "ddsim",
        "--compactFile", compact,
        "--outputFile", out,
        "-N", str(events),
        "--enableGun",
        "--gun.particle", "mu-",
        "--gun.energy", "10*GeV",
        "--gun.position", f"{x}*mm {y}*mm -100*mm",
        "--gun.direction", "0 0 1",
        "--part.userParticleHandler=",
        "--physicsList", "QGSP_BERT",
    ]
    log = out + ".log"
    with open(log, "w") as fh:
        rc = subprocess.call(cmd, stdout=fh, stderr=subprocess.STDOUT)
    if rc != 0:
        sys.exit(f"ddsim failed (rc={rc}), see {log}")


_DECODER = None


def decoder(compact: str):
    """CellID decoder for the compact file. Detector is a singleton, so the
    geometry can only be loaded once per process."""
    global _DECODER
    if _DECODER is None:
        import dd4hep

        det = dd4hep.Detector.getInstance()
        det.fromXML(compact)
        _DECODER = det.readout(READOUT).idSpec().decoder()
    return _DECODER


def read_hits(compact: str, path: str, nlayers: int):
    """Return (efficiency, Counter of (x, y) cells).

    Efficiency is the fraction of (event, layer) pairs that produced at least one
    hit -- i.e. 'did the muon leave a signal in this plane at all'. Counting raw
    hits instead would overshoot 100%, because delta rays sprinkle extra hits in
    neighbouring pads.
    """
    from podio.root_io import Reader

    dec = decoder(compact)
    cells, fired, events = collections.Counter(), 0, 0
    for event in Reader(path).get("events"):
        events += 1
        layers = set()
        for hit in event.get(READOUT):
            cid = hit.getCellID()
            cells[(dec.get(cid, "x"), dec.get(cid, "y"))] += 1
            layers.add(dec.get(cid, "layer"))
        fired += len(layers)
    return fired / (events * nlayers), cells


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--compact", default=DEFAULT_COMPACT)
    ap.add_argument("--events", type=int, default=20)
    ap.add_argument("--layers", type=int, default=15, help="expected hits per event")
    ap.add_argument("--x", type=float, help="single target x [mm]")
    ap.add_argument("--y", type=float, help="single target y [mm]")
    args = ap.parse_args()

    if args.x is not None and args.y is not None:
        targets = [(args.x, args.y, "requested point")]
    else:
        targets = DEFAULT_TARGETS

    workdir = tempfile.mkdtemp(prefix="dead_zone_")
    print(f"compact : {args.compact}")
    print(f"muons   : {args.events} x 10 GeV through {args.layers} layers")
    print(f"workdir : {workdir}\n")
    print(f"{'target':<52}{'plane eff':>10}   dominant cells")

    for x, y, label in targets:
        out = os.path.join(workdir, f"mu_{x}_{y}.root".replace("-", "m"))
        simulate(args.compact, x, y, args.events, out)
        eff, cells = read_hits(args.compact, out, args.layers)
        top = ", ".join(f"{c}:{n}" for c, n in cells.most_common(2))
        print(f"{label:<52}{eff:>9.0%}   {top}")

    print("\nA target on silicon must read 100%. A target in a dead region must read\n"
          "well below it; the residual efficiency is muons scattering into a\n"
          "neighbouring pad deeper in the stack.")


if __name__ == "__main__":
    main()

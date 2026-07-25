#!/usr/bin/env python3
"""Rewrite the pad maps onto the pad grid of the simulated geometry.

The maps are the ``(chip, channel) -> (x, y)`` tables the test-beam analysis and
ChannelMapper use. They were written on the nominal design grid: a round 5.5 mm
pitch, 7.6 mm across the boundary between two sensors. The geometry now uses the
measured sensor -- a 5.52 mm diode on a 5.53 mm pixel pitch, with a 0.61 mm
inactive rim and no guard ring -- which puts the pads up to 0.43 mm away from
where the nominal table has them. Small enough that ChannelMapper's nearest-pad
matching never noticed, but it means the two tables disagree, and only one of
them can be the geometry the analysis assumes.

This script maps the old columns onto the new ones and rewrites the files in
place. Only coordinates change: the chip/channel assignment, the row order and
the header are copied through untouched, so anything keyed on (chip, channel) --
the masking lists, the MIP calibration -- is unaffected.

The transformation is a linear interpolation through the 32 old columns and the
32 new ones. On a pad centre it is exact by construction; the ``x0``/``y0``
reference columns, which are not pad centres, are carried along by the same
interpolation, so a chip's reference stays at the same place relative to its
pads. Both consumers of the file (ChannelMapper, event_viewer._pad_map) read
those two columns and discard them, but they are kept consistent anyway.

The new grid is read from the compact file, never hardcoded, so re-running this
after a geometry change is all it takes.

Usage:
    source /cvmfs/sw.hsf.org/key4hep/setup.sh -r 2026-04-08
    export LD_LIBRARY_PATH=$PWD/install/lib64:$PWD/install/lib:$LD_LIBRARY_PATH
    python mappings/regenerate_pad_maps.py            # rewrite in place
    python mappings/regenerate_pad_maps.py --dry-run  # just report the shift
"""

from __future__ import annotations

import argparse
import os
import shutil
import sys

import numpy as np

REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
COMPACT = os.path.join(REPO_ROOT, "simulation", "geometry", "SND_compact.xml")
MAPPINGS = os.path.join(REPO_ROOT, "mappings")

MAP_FILES = [
    "fev10_rotate_chip_channel_x_y_mapping.txt",
    "fev11_cob_good_rotate_chip_channel_x_y_mapping.txt",
]

CM_TO_MM = 10.0  # DD4hep works in cm, the maps are written in mm


def geometry_columns(compact: str) -> np.ndarray:
    """Pad centres of one layer along x [mm], from the compact file.

    A pad's centre is its position inside its sensor's pad array, plus the
    position of that sensor: the array is centred on the sensor, so the sensor's
    inactive rim shifts nothing, it only sets how far apart two sensors are.
    """
    import dd4hep

    det = dd4hep.Detector.getInstance()
    if det.state() != dd4hep.Detector.READY:
        det.fromXML(compact)
    c = det.constantAsDouble

    pitch = c("Ecal_CellSizeX")
    active = c("Ecal_WaferActiveX")
    n_pads = int(c("Ecal_NPadsPerWaferX"))
    sensor = c("Ecal_WaferSizeX") + c("Ecal_WaferGapX")   # sensor to sensor
    n_sensors = int(c("Ecal_NWafersX")) * int(c("Ecal_NASUsX"))

    columns = []
    for s in range(n_sensors):
        centre = (s - (n_sensors - 1) / 2.0) * sensor
        for p in range(n_pads):
            columns.append(centre - active / 2.0 + (p + 0.5) * pitch)
    return np.asarray(sorted(columns)) * CM_TO_MM


def read_rows(path: str):
    """``(header, [(chip, x0, y0, channel, x, y)])`` of a pad map file."""
    header, rows = None, []
    for n, line in enumerate(open(path)):
        parts = line.split()
        if not parts:
            continue
        try:
            chip = int(parts[0])
        except ValueError:
            if header is None:
                header = line.rstrip("\n")
            continue
        rows.append((chip, float(parts[1]), float(parts[2]),
                     int(parts[3]), float(parts[4]), float(parts[5])))
    return header, rows


def fmt(value: float) -> str:
    """Match the style of the existing files: plain decimals, no trailing zeros."""
    text = f"{value:.4f}".rstrip("0").rstrip(".")
    return "0" if text in ("", "-0") else text


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--compact", default=COMPACT)
    ap.add_argument("--dry-run", action="store_true",
                    help="report the shift without touching the files")
    ap.add_argument("--no-backup", action="store_true",
                    help="do not leave a .nominal copy of the old file")
    args = ap.parse_args()

    new_cols = geometry_columns(args.compact)
    print(f"geometry: {new_cols.size} columns, "
          f"{new_cols[0]:.4f} .. {new_cols[-1]:.4f} mm")

    for name in MAP_FILES:
        path = os.path.join(MAPPINGS, name)
        header, rows = read_rows(path)
        old_cols = np.asarray(sorted({round(r[4], 4) for r in rows}
                                     | {round(r[5], 4) for r in rows}))
        if old_cols.size != new_cols.size:
            sys.exit(f"{name}: map has {old_cols.size} columns, "
                     f"geometry has {new_cols.size}")

        def remap(v: float) -> float:
            return float(np.interp(v, old_cols, new_cols))

        shift = max(abs(remap(c) - c) for c in old_cols)
        print(f"{name}: {len(rows)} pads, largest move {shift:.4f} mm")
        if args.dry_run:
            continue

        if not args.no_backup:
            shutil.copyfile(path, path + ".nominal")
        with open(path, "w") as out:
            if header:
                out.write(header + "\n")
            for chip, x0, y0, channel, x, y in rows:
                out.write(f"{chip} {fmt(remap(x0))} {fmt(remap(y0))} "
                          f"{channel} {fmt(remap(x))} {fmt(remap(y))}\n")

    if args.dry_run:
        print("\ndry run, nothing written")


if __name__ == "__main__":
    main()

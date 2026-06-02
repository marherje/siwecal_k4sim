#!/usr/bin/env python3
"""Diagnostic: truth muon trajectory (PDG-filtered SimCalorimeterHit positions)
vs ACTS reco trajectory, spanning the full SND detector.

Reads per-detector hits from timewindows.edm4hep.root, tags each hit by the
PDG of its highest-energy contribution (from the parallel '{det}ContribPDGs'
frame parameter), and overlays them with ACTS track states from ShipHits.root.
Produces per-window full_detector_w*.pdf and a quantitative MTC curvature table.

Note: DDG4 in key4hep-2026-02-01 does not transfer CaloHitContribution.stepPosition,
so we read hit.getPosition() (per-hit, written by the sensitive actions) instead.
All positions are in mm.  ACTSTrackStates x, y, z also in mm.
"""

import argparse
import math
import os
import sys
from pathlib import Path
from collections import defaultdict

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import ROOT

HERE = Path(__file__).resolve().parent
PROJECT = HERE.parent.parent
sys.path.insert(0, str(PROJECT / "simulation" / "geometry"))
from parse_geometry import SNDGeometry  # type: ignore  # noqa: E402

# ---------------------------------------------------------------------------
# Geometry constants
# ---------------------------------------------------------------------------
GEO = SNDGeometry()
MTC_Z_CENTERS = GEO.mtc_station_z_centers          # [mm]
MTC_TOTAL_LEN = GEO.mtc_n_layers * GEO.mtc_layer_thick
MTC_Z_MIN = min(MTC_Z_CENTERS) - MTC_TOTAL_LEN / 2.0 - 5.0
MTC_Z_MAX = max(MTC_Z_CENTERS) + MTC_TOTAL_LEN / 2.0 + 5.0
TAN5 = math.tan(5.0 * math.pi / 180.0)  # SciFi stereo angle tangent

# Detector display colours for shading bands
_DET_COLORS = {
    "SiTarget": "steelblue",
    "SiPad":    "darkorange",
    "MTCSciFi": "forestgreen",
    "MTCScint": "orchid",
}
# Muon scatter colours per detector (|PDG|==13 hits)
_MUON_COLORS = {
    "SiTarget": "cyan",
    "SiPad":    "orange",
    "MTCSciFi": "lime",
    "MTCScint": "magenta",
}

# cellID layouts (from simulation/geometry/SND_compact.xml):
#   SiTarget:  system:8,layer:8,slice:4,plane:1,column:2,row:2,strip:14
#              -> plane bit at offset 20 (1 bit)
#              plane=0 -> X-strip (measures x), plane=1 -> Y-strip (measures y)
#   MTC:       system:8,station:2,layer:8,slice:4,plane:2,strip:14,x:9,y:9
#              -> plane field at offset 22 (2 bits)
#              plane=0 -> SciFi U (+5 deg stereo), plane=1 -> SciFi V (-5 deg), plane=2 -> Scint
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

# ---------------------------------------------------------------------------
# Per-detector contribution reader (podio)
# ---------------------------------------------------------------------------
def _open_podio_reader(edm4hep_file):
    try:
        from podio import reading
        return reading.get_reader(edm4hep_file)
    except Exception:
        import podio
        return podio.root_io.Reader(edm4hep_file)


def read_contribs_by_pdg(edm4hep_file, det_name, mock_pdgs=False):
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
    out = {}
    warned = False
    for iframe, frame in enumerate(reader.get("events")):
        try:
            hits = frame.get(f"{det_name}HitsWindowed")
            contribs = frame.get(f"{det_name}HitsWindowed_Contributions")
        except Exception:
            continue
        if hits is None or contribs is None:
            continue

        # PDG list from frame parameter.  podio returns None when absent, an int
        # for single-element vectors, and a list otherwise -- normalise to a list.
        if mock_pdgs:
            pdg_list = [13] * len(contribs)
        else:
            raw = frame.get_parameter(f"{det_name}ContribPDGs")
            if raw is None:
                if not warned:
                    print(f"[warn] {det_name}ContribPDGs parameter absent "
                          f"-- tagging all hits as PDG=0 (old file?)")
                    warned = True
                pdg_list = [0] * len(contribs)
            elif isinstance(raw, int):
                pdg_list = [raw]
            else:
                pdg_list = list(raw)

        # Build ObjectID-index -> PDG lookup so we can find each linked
        # contribution's PDG via the relation on the hit.
        oid_to_pdg = {}
        for ci, c in enumerate(contribs):
            if ci < len(pdg_list):
                oid_to_pdg[c.getObjectID().index] = int(pdg_list[ci])

        by_pdg = defaultdict(list)
        for hit in hits:
            # Highest-energy contribution wins for the per-hit tag.
            best_pdg, best_E = 0, -1.0
            for c in hit.getContributions():
                p = oid_to_pdg.get(c.getObjectID().index, 0)
                e = float(c.getEnergy())
                if e > best_E:
                    best_E, best_pdg = e, p
            pos = hit.getPosition()  # mm
            # Approximate per-hit timing/energy: take the first contribution's time
            # and the hit's total energy.  Time is window-local (set by splitter).
            try:
                t = float(next(iter(hit.getContributions())).getTime())
            except StopIteration:
                t = 0.0
            E = float(hit.getEnergy())
            plane = _hit_plane(int(hit.getCellID()), det_name)
            by_pdg[best_pdg].append(
                (float(pos.x), float(pos.y), float(pos.z), t, E, plane)
            )

        out[iframe] = dict(by_pdg)
    return out


# ---------------------------------------------------------------------------
# Per-panel hit selection (per detector) + SciFi U/V stereo recovery
# ---------------------------------------------------------------------------
# Each `pts` argument below is a list of (x, y, z, t, E, plane) tuples.

def hits_for_xpanel(pts, det_name):
    """Return [(x_meaningful, z, E), ...] suitable for the z-x panel.

    SiTarget: only X-strip hits (plane=0); position.x is the muon x measurement.
    MTC SciFi: U+V stereo recovery via recover_scifi_x.
    SiPad / MTCScint: pad/grid centre — position.x is the muon x.
    """
    if det_name == "SiTarget":
        return [(x, z, E) for (x, y, z, t, E, plane) in pts if plane == 0]
    if det_name == "MTCSciFi":
        return recover_scifi_x(pts)
    return [(x, z, E) for (x, y, z, t, E, plane) in pts]


def hits_for_ypanel(pts, det_name):
    """Return [(y_meaningful, z, E), ...] suitable for the z-y panel.

    SiTarget: only Y-strip hits (plane=1).
    MTC SciFi: position.y is the energy-weighted step Y (real muon y, written
               by SND_SciFiAction) on every plane, so all SciFi hits are valid.
    SiPad / MTCScint: pad/grid centre.
    """
    if det_name == "SiTarget":
        return [(y, z, E) for (x, y, z, t, E, plane) in pts if plane == 1]
    return [(y, z, E) for (x, y, z, t, E, plane) in pts]


def recover_scifi_x(pts, layer_bin_mm=5.0):
    """Recover muon x at each MTC SciFi layer by averaging U (plane=0) and V
    (plane=1) strip centres at the same z.  Unpaired hits (only U or only V in
    a layer) are linearly interpolated in z from the surrounding paired layers;
    if no paired layer exists on both sides of an unpaired hit, the hit is
    excluded.  Returns [(gx, gz, E), ...] sorted by z.
    """
    if not pts:
        return []
    layer = lambda z: int(round(z / layer_bin_mm))
    u_by_bin = defaultdict(list)
    v_by_bin = defaultdict(list)
    for x, y, z, t, E, plane in pts:
        if plane == 0:
            u_by_bin[layer(z)].append((x, z, E))
        elif plane == 1:
            v_by_bin[layer(z)].append((x, z, E))

    paired = []          # (gx, gz, E_total)
    unpaired_singles = []  # (x_strip, z, E)
    for bin_id in sorted(set(u_by_bin.keys()) | set(v_by_bin.keys())):
        u_grp = u_by_bin.get(bin_id, [])
        v_grp = v_by_bin.get(bin_id, [])
        if u_grp and v_grp:
            u_E = sum(h[2] for h in u_grp) or 1e-30
            v_E = sum(h[2] for h in v_grp) or 1e-30
            u_x = sum(h[0] * h[2] for h in u_grp) / u_E
            v_x = sum(h[0] * h[2] for h in v_grp) / v_E
            u_z = sum(h[1] * h[2] for h in u_grp) / u_E
            v_z = sum(h[1] * h[2] for h in v_grp) / v_E
            paired.append(((u_x + v_x) / 2.0, (u_z + v_z) / 2.0, u_E + v_E))
        else:
            for h in (u_grp or v_grp):
                unpaired_singles.append(h)

    out = list(paired)
    if paired and unpaired_singles:
        paired.sort(key=lambda p: p[1])
        zs = np.array([p[1] for p in paired])
        xs = np.array([p[0] for p in paired])
        for _, z, E in unpaired_singles:
            # np.interp clamps at the boundaries, so reject any hit outside the
            # paired range — we only want true interpolation, never extrapolation.
            if z < zs[0] or z > zs[-1]:
                continue
            out.append((float(np.interp(z, zs, xs)), float(z), E))

    return sorted(out, key=lambda p: p[1])


def _track_muon_path(coord_z_e, z_bin_mm=5.0):
    """[(coord, z, E), ...] -> [(coord, z), ...] tracking the muon trajectory.

    Groups hits into ~5 mm z-bins (≈ one MTC SciFi layer pair) and picks ONE
    representative point per bin:

      * one candidate in a bin  -> use it as is
      * first multi-candidate bin            -> highest-energy hit
      * second multi-candidate bin           -> nearest in coord to the previous point
      * any later multi-candidate bin        -> nearest in coord to the linear
        extrapolation of the last two path points to that bin's z

    This is needed because a single SciFi plane can have multiple strips lit
    in the same layer when delta-rays / brem-induced electrons accompany the
    muon.  ddsim does not save those secondaries as separate MCParticles, so
    they back-link to the primary muon and the PDG=13 filter cannot separate
    them.  Tracking along the trajectory picks the correct strip per layer.
    """
    if not coord_z_e:
        return []
    bins = defaultdict(list)
    for c, z, E in coord_z_e:
        bins[int(round(z / z_bin_mm))].append((c, z, E))

    path = []
    for li in sorted(bins.keys()):
        cands = bins[li]
        if len(cands) == 1:
            c, z, _ = cands[0]
        elif not path:
            c, z, _ = max(cands, key=lambda h: h[2])
        else:
            avg_z = sum(h[1] for h in cands) / len(cands)
            if len(path) == 1:
                pred = path[-1][0]
            else:
                (c0, z0), (c1, z1) = path[-2], path[-1]
                pred = (c1 + (avg_z - z1) * (c1 - c0) / (z1 - z0)
                        if z1 != z0 else c1)
            c, z, _ = min(cands, key=lambda h: abs(h[0] - pred))
        path.append((c, z))
    return path


def extract_muon_paths(by_pdg, det_name):
    """Return {'zx': [(x, z), ...], 'zy': [(y, z), ...]} for muon hits only.

    Per-layer selection via nearest-neighbour tracking (_track_muon_path)
    rather than energy-weighted averaging, because energy-weighting is fooled
    by delta-ray strips that often deposit more energy than the muon's own
    strip in the same layer.
    """
    pts = by_pdg.get(13, []) + by_pdg.get(-13, [])
    if not pts:
        return {'zx': [], 'zy': []}
    return {
        'zx': _track_muon_path(hits_for_xpanel(pts, det_name)),
        'zy': _track_muon_path(hits_for_ypanel(pts, det_name)),
    }


# ---------------------------------------------------------------------------
# Reco track states reader (RNTuple)
# ---------------------------------------------------------------------------
def read_reco_track_states(ship_hits_file):
    """Returns dict {window_id: {track_id: [(x_mm, y_mm, z_mm, loc0, tilt, chi2), ...]}}.

    chi2 here is the PER-STATE innovation chi² (≈ calibratedSize for a good fit:
    ~1 for 1D strip surfaces, ~2 for 2D pad/combined surfaces). For the TRACK-level
    chi²/ndf goodness-of-fit summary, use read_reco_tracks_meta().

    Uses RDF.FromRNTuple → AsNumpy to avoid RNTuple teardown segfaults.
    """
    out = defaultdict(lambda: defaultdict(list))
    try:
        df = ROOT.RDF.FromRNTuple("ACTSTrackStates", ship_hits_file)
        cols = ["window_id", "track_id", "x", "y", "z", "loc0", "tilt", "chi2"]
        data = df.AsNumpy(cols)
    except Exception as e:
        print(f"[reco] cannot open ACTSTrackStates: {e}")
        return out
    n = len(data["window_id"])
    for i in range(n):
        wid = int(data["window_id"][i])
        tid = int(data["track_id"][i])
        out[wid][tid].append((
            float(data["x"][i]),
            float(data["y"][i]),
            float(data["z"][i]),
            float(data["loc0"][i]),
            float(data["tilt"][i]),
            float(data["chi2"][i]),
        ))
    return out


def read_reco_tracks_meta(ship_hits_file, track_collection="ACTSTracks"):
    """Returns dict {(window_id, track_id): {'chi2': float, 'ndf': int, 'chi2_per_ndf': float}}.

    Source: the per-track RNTuple written by EDM4HEP2RNTuple, where
        ndf  = Σ calibratedSize − n_fit_params (5 helix params)
        chi² = Σ per-state innovation chi² (ACTS filter pass)
    so chi²/ndf is the standard goodness-of-fit ratio (≈1 = good).
    """
    out = {}
    try:
        df = ROOT.RDF.FromRNTuple(track_collection, ship_hits_file)
        data = df.AsNumpy(["window_id", "track_id", "chi2", "ndf"])
    except Exception as e:
        print(f"[reco] cannot open {track_collection}: {e}")
        return out
    n = len(data["window_id"])
    for i in range(n):
        wid  = int(data["window_id"][i])
        tid  = int(data["track_id"][i])
        chi2 = float(data["chi2"][i])
        ndf  = int(data["ndf"][i])
        out[(wid, tid)] = {
            "chi2":         chi2,
            "ndf":          ndf,
            "chi2_per_ndf": chi2 / ndf if ndf > 0 else float("nan"),
        }
    return out


def stereo_pair_track_states(states):
    """Compute (x, y, z) from loc0_U / loc0_V stereo pairs instead of ACTS smoothed values.

    U/V SciFi planes in one MTC layer are ~2.35 mm apart in z. For each
    adjacent (U, V) pair within a 5 mm z-window:
        gx = (loc0_U + loc0_V) / 2
        gy = (loc0_V - loc0_U) / (2 * tan5°)
    Non-stereo states (|tilt| ≈ 0, i.e. scintillator) fall back to smoothed x/y.
    Returns [(x_mm, y_mm, z_mm), ...] sorted by z.
    """
    states_sorted = sorted(states, key=lambda s: s[2])
    pts = []
    i = 0
    print("CHECK")
    while i < len(states_sorted):
        x, y, z, loc0, tilt = states_sorted[i]
        if abs(tilt) > 0.01 and i + 1 < len(states_sorted):
            _, _, z2, loc02, tilt2 = states_sorted[i + 1]
            print("check condition...")
            if abs(z - z2) < 5.0 and tilt * tilt2 < 0:
                loc0_U = loc0  if tilt > 0 else loc02
                loc0_V = loc02 if tilt > 0 else loc0
                gx = (loc0_U + loc0_V) / 2.0
                gy = (loc0_V - loc0_U) / (2.0 * TAN5)
                gz = (z + z2) / 2.0
                pts.append((gx, gy, gz))
                i += 2
                print("PAIRED: ", gx, gy, gz)
                continue
        pts.append((x, y, z))
        i += 1
    pts.sort(key=lambda p: p[2])
    return pts


def raw_track_states_xyz(states, use_loc0_y=False):
    """Return (x_mm, y_mm, z_mm) sorted by z.

    use_loc0_y=True: recover x/y from SciFi loc0_U/loc0_V stereo pairing.
    use_loc0_y=False (default): use ACTS smoothed x/y directly.
    """
    if use_loc0_y:
        return stereo_pair_track_states(states)
    return [(x, y, z) for (x, y, z, _, _) in sorted(states, key=lambda s: s[2])]


# ---------------------------------------------------------------------------
# Full-detector 2-panel plot
# ---------------------------------------------------------------------------
def plot_full_detector(window_id, contribs_by_det, muon_paths_by_det,
                       reco_tracks, reco_meta, outdir, use_loc0_y=False):
    """2-panel (z-x, z-y) plot spanning all detectors for one window.

    Per-panel hit selection (so SiTarget X/Y-strip alternation no longer
    zigzags, and MTC SciFi z-x uses U+V stereo-recovered x):
      z-x panel: hits_for_xpanel(...) per detector
      z-y panel: hits_for_ypanel(...) per detector

    contribs_by_det : {det: {pdg: [(x, y, z, t, E, plane), ...]}}
    muon_paths_by_det: {det: {'zx': [(x, z)], 'zy': [(y, z)]}}
    reco_tracks      : {track_id: [(x, y, z, loc0, tilt, chi2), ...]}
    reco_meta        : {(window_id, track_id): {'chi2', 'ndf', 'chi2_per_ndf'}}
    """
    fig, (ax_zx, ax_zy) = plt.subplots(2, 1, figsize=(14, 10), sharex=True)

    # Derive sub-detector z extents from data (data-driven, no hardcodes)
    det_z_extents = {}
    for det, by_pdg in contribs_by_det.items():
        all_z = [z for pts in by_pdg.values() for (x, y, z, t, E, p) in pts]
        if all_z:
            det_z_extents[det] = (min(all_z), max(all_z))

    # Shade sub-detector bands (behind everything)
    for det, (zlo, zhi) in det_z_extents.items():
        color = _DET_COLORS.get(det, "gray")
        for ax in (ax_zx, ax_zy):
            ax.axvspan(zlo, zhi, alpha=0.06, color=color, label=f"_{det}_band")
            ax.text((zlo + zhi) / 2.0, 1.01, det, transform=ax.get_xaxis_transform(),
                    ha="center", va="bottom", fontsize=7, color=color)

    # MTC station fine-grain shading
    for zc in MTC_Z_CENTERS:
        for ax in (ax_zx, ax_zy):
            ax.axvspan(zc - MTC_TOTAL_LEN / 2.0, zc + MTC_TOTAL_LEN / 2.0,
                       alpha=0.05, color="gray")

    # Per-detector scatter: panel-correct hit selection.
    for det, by_pdg in contribs_by_det.items():
        muon_color = _MUON_COLORS.get(det, "yellow")
        pts_all = [pt for pts in by_pdg.values() for pt in pts]
        pts_mu  = by_pdg.get(13, []) + by_pdg.get(-13, [])

        zx_all = hits_for_xpanel(pts_all, det)
        zy_all = hits_for_ypanel(pts_all, det)
        zx_mu  = hits_for_xpanel(pts_mu, det)
        zy_mu  = hits_for_ypanel(pts_mu, det)

        if zx_all:
            ax_zx.scatter([p[1] for p in zx_all], [p[0] for p in zx_all],
                          s=2, c="lightgray", alpha=0.25, marker=".",
                          rasterized=True)
        if zy_all:
            ax_zy.scatter([p[1] for p in zy_all], [p[0] for p in zy_all],
                          s=2, c="lightgray", alpha=0.25, marker=".",
                          rasterized=True)
        if zx_mu:
            ax_zx.scatter([p[1] for p in zx_mu], [p[0] for p in zx_mu],
                          s=8, c=muon_color, alpha=0.7, marker="o", zorder=3,
                          label=f"{det} μ x ({len(zx_mu)})")
        if zy_mu:
            ax_zy.scatter([p[1] for p in zy_mu], [p[0] for p in zy_mu],
                          s=8, c=muon_color, alpha=0.7, marker="o", zorder=3,
                          label=f"{det} μ y ({len(zy_mu)})")

    # Unified muon truth path per panel: concat detectors in _TRUTH_PATH_DETS,
    # sort by z.  MTCScint is excluded (coarse pad granularity).
    unified_zx = sorted(
        [pt for det in _TRUTH_PATH_DETS if det in muon_paths_by_det
         for pt in muon_paths_by_det[det]['zx']],
        key=lambda p: p[1])
    unified_zy = sorted(
        [pt for det in _TRUTH_PATH_DETS if det in muon_paths_by_det
         for pt in muon_paths_by_det[det]['zy']],
        key=lambda p: p[1])
    if unified_zx:
        ax_zx.plot([p[1] for p in unified_zx], [p[0] for p in unified_zx],
                   "k-o", ms=4, lw=1.5, zorder=5,
                   label=f"Truth μ x ({len(unified_zx)} bins)")
    if unified_zy:
        ax_zy.plot([p[1] for p in unified_zy], [p[0] for p in unified_zy],
                   "k-o", ms=4, lw=1.5, zorder=5,
                   label=f"Truth μ y ({len(unified_zy)} bins)")

    # ACTS reco track states
    reco_colors = ["red", "blue", "green", "magenta"]
    for it, (tid, raw_states) in enumerate(reco_tracks.items()):
        pts = raw_track_states_xyz([x[:-1] for x in raw_states], use_loc0_y=use_loc0_y)
        if not pts:
            continue
        rz = [p[2] for p in pts]
        rx = [p[0] for p in pts]
        ry = [p[1] for p in pts]
        c = reco_colors[it % len(reco_colors)]
        meta = reco_meta.get((window_id, tid), {})
        chi2_ndf = meta.get("chi2_per_ndf", float("nan"))
        ax_zx.plot(rz, rx, "-o", ms=3, c=c, lw=1.2, zorder=4,
                   label=f"Reco track {tid}: chi2/ndf={chi2_ndf:.2f}")
        ax_zy.plot(rz, ry, "-o", ms=3, c=c, lw=1.2, zorder=4)

    # Axis limits from data — fall back to MTC range if no hits
    all_zs_plot = [z for by_pdg in contribs_by_det.values()
                   for pts in by_pdg.values() for (x, y, z, t, E, p) in pts]
    all_zs_plot += [p[1] for p in unified_zx] + [p[1] for p in unified_zy]
    if all_zs_plot:
        zlo_plot = min(all_zs_plot) - 20.0
        zhi_plot = max(all_zs_plot) + 20.0
    else:
        zlo_plot, zhi_plot = MTC_Z_MIN, MTC_Z_MAX

    for ax, ylabel in [(ax_zx, "x [mm]"), (ax_zy, "y [mm]")]:
        ax.set_xlim(zlo_plot, zhi_plot)
        ax.axhline(0, lw=0.4, color="gray")
        ax.set_ylabel(ylabel)
        ax.grid(alpha=0.3)
        ax.legend(loc="best", fontsize=7, ncol=2)

    ax_zy.set_xlabel("z (beam) [mm]")
    ax_zx.set_title(f"Window {window_id}: full detector z–x")
    ax_zy.set_title(f"Window {window_id}: full detector z–y")

    outdir.mkdir(parents=True, exist_ok=True)
    outfile = outdir / f"full_detector_w{window_id:03d}.pdf"
    fig.tight_layout()
    fig.savefig(outfile)
    plt.close(fig)
    return outfile


# ---------------------------------------------------------------------------
# Quantitative summary (MTC curvature)
# ---------------------------------------------------------------------------
def fit_curvature_in_iron(zs, xs):
    """Parabolic fit x(z) ~ x0 + slope*Δz + 0.5*curv*Δz².
    Returns dict with radius_mm = 1/curv, or None if < 3 points.
    """
    if len(zs) < 3:
        return None
    z = np.asarray(zs, dtype=float)
    x = np.asarray(xs, dtype=float)
    z0 = z.mean()
    c2, c1, c0 = np.polyfit(z - z0, x, 2)
    if abs(c2) < 1e-20:
        return None
    R = 1.0 / (2.0 * c2)
    return {"x0": float(c0), "slope": float(c1), "curv": float(2.0 * c2),
            "radius_mm": float(R)}


def momentum_from_radius(radius_mm, B_T=1.7):
    """p [GeV] = 0.3 * B[T] * r[m].  r in mm → /1000."""
    return 0.3 * B_T * (radius_mm / 1000.0)


def summarize(window_id, truth_zx_path, reco_tracks, reco_meta, use_loc0_y=False):
    """Compute MTC parabolic curvature deviation between truth and reco.

    truth_zx_path: list of (x, z) — the unified MTC SciFi U/V-recovered
    muon x-vs-z trajectory (plus SiPad / SiTarget x-strip points where present).
    """
    if not truth_zx_path:
        return None
    tz = np.array([p[1] for p in truth_zx_path])
    tx = np.array([p[0] for p in truth_zx_path])
    mtc_mask = (tz >= MTC_Z_MIN) & (tz <= MTC_Z_MAX)
    if mtc_mask.sum() < 3:
        return None
    tz_mtc = tz[mtc_mask]
    tx_mtc = tx[mtc_mask]
    truth_fit = fit_curvature_in_iron(tz_mtc.tolist(), tx_mtc.tolist())

    results = []
    for tid, raw_states in reco_tracks.items():
        pts = raw_track_states_xyz([x[:-1] for x in raw_states], use_loc0_y=use_loc0_y)
        if not pts:
            continue
        rz = np.array([p[2] for p in pts])
        rx = np.array([p[0] for p in pts])
        in_mtc = (rz >= MTC_Z_MIN) & (rz <= MTC_Z_MAX)
        if in_mtc.sum() < 3:
            continue
        rz_m = rz[in_mtc]
        rx_m = rx[in_mtc]
        tx_interp = np.interp(rz_m, tz_mtc, tx_mtc)
        dx = rx_m - tx_interp
        reco_fit = fit_curvature_in_iron(rz_m.tolist(), rx_m.tolist())
        meta = reco_meta.get((window_id, tid), {})
        results.append({
            "track_id":      tid,
            "n_in_mtc":      int(in_mtc.sum()),
            "x_dev_mean":    float(np.mean(np.abs(dx))),
            "x_dev_max":     float(np.max(np.abs(dx))),
            "x_first_mtc":   float(rx_m[0]),
            "x_last_mtc":    float(rx_m[-1]),
            "x_first_truth": float(tx_interp[0]),
            "x_last_truth":  float(tx_interp[-1]),
            "reco_fit":      reco_fit,
            "truth_fit":     truth_fit,
            "chi2":          meta.get("chi2",         float("nan")),
            "ndf":           meta.get("ndf",          0),
            "chi2/ndf":      meta.get("chi2_per_ndf", float("nan")),
        })
    return {"window_id": window_id,
            "n_truth_mtc": int(mtc_mask.sum()),
            "tracks": results}


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--truth", default=str(HERE / "timewindows.edm4hep.root"),
                    help="EDM4HEP file with windowed hits (timewindows.edm4hep.root)")
    ap.add_argument("--reco",  default=str(HERE / "ShipHits.root"),
                    help="ShipHits.root RNTuple with ACTSTrackStates")
    ap.add_argument("--max-windows", type=int, default=10,
                    help="Plot only the first N windows")
    ap.add_argument("--outdir", default=str(HERE / "plots_curvature"),
                    help="Output directory for plots")
    ap.add_argument("--mock-pdgs", action="store_true",
                    help="DEV-only: label every contribution PDG=13 to smoke-test "
                         "plotting before the C++ PDG-parameter pipeline is merged.")
    ap.add_argument("--use-loc0-y", action="store_true",
                    help="Recover reco x/y from SciFi loc0_U/loc0_V stereo pairing "
                         "instead of using ACTS smoothed x/y directly.")
    args = ap.parse_args()

    print(f"[diag] truth:       {args.truth}")
    print(f"[diag] reco:        {args.reco}")
    print(f"[diag] mock-pdgs:   {args.mock_pdgs}")
    print(f"[diag] use-loc0-y:  {args.use_loc0_y}")
    if args.mock_pdgs:
        print("[diag] WARNING: --mock-pdgs active — PDG labels are fabricated!")

    DET_NAMES = ("SiTarget", "SiPad", "MTCSciFi", "MTCScint")

    # Read all detectors up-front (one pass per detector through the file)
    print("[diag] reading hit positions tagged by per-contribution PDG …")
    all_contribs = {}
    for det in DET_NAMES:
        all_contribs[det] = read_contribs_by_pdg(args.truth, det,
                                                   mock_pdgs=args.mock_pdgs)
        n_frames = len(all_contribs[det])
        n_total  = sum(sum(len(v) for v in bw.values())
                       for bw in all_contribs[det].values())
        print(f"[diag]   {det}: {n_frames} frames, {n_total} hits")

    reco      = read_reco_track_states(args.reco)
    reco_meta = read_reco_tracks_meta(args.reco)
    print(f"[diag] reco windows: {len(reco)}  meta tracks: {len(reco_meta)}")

    # Union of all window IDs that have any data
    all_windows_set = set()
    for det in DET_NAMES:
        all_windows_set |= set(all_contribs[det].keys())
    all_windows_set |= set(reco.keys())
    all_windows = sorted(all_windows_set)[: args.max_windows]

    outdir = Path(args.outdir)
    summaries = []

    for w in all_windows:
        # Per-window contribs dict: {det: {pdg: [(x,y,z,t,E,plane), ...]}}
        contribs_w = {det: all_contribs[det].get(w, {}) for det in DET_NAMES}

        # Per-detector muon truth paths split into {'zx': [(x,z)], 'zy': [(y,z)]}
        muon_paths_w = {det: extract_muon_paths(contribs_w[det], det)
                        for det in DET_NAMES}

        # Unified z-x truth (used for the MTC curvature fit) — concat across
        # _TRUTH_PATH_DETS only (excludes MTCScint).
        unified_zx = sorted(
            [pt for det in _TRUTH_PATH_DETS
             for pt in muon_paths_w[det]['zx']],
            key=lambda p: p[1])

        reco_trks = reco.get(w, {})
        outfile = plot_full_detector(w, contribs_w, muon_paths_w, reco_trks,
                                     reco_meta, outdir,
                                     use_loc0_y=args.use_loc0_y)
        zx_counts = ", ".join(f"{det}={len(muon_paths_w[det]['zx'])}"
                              for det in DET_NAMES)
        print(f"[diag] window {w}: unified_zx_truth={len(unified_zx)} bins "
              f"(zx mu bins per det: {zx_counts}), "
              f"reco_tracks={len(reco_trks)} → {outfile.name}")

        s = summarize(w, unified_zx, reco_trks, reco_meta,
                      use_loc0_y=args.use_loc0_y)
        if s and s["tracks"]:
            summaries.append(s)

    print("\n===== MTC PARABOLIC-FIT CURVATURE COMPARISON =====")
    print("Reco x(z) and truth-muon x(z) inside MTC are each fit to a parabola")
    print("x = c0 + c1*Δz + 0.5*curv*Δz².  Effective curvature radius R = 1/curv [mm].")
    print("Equivalent average momentum p_eff [GeV] = 0.3 * B[T=1.7] * R[m] (only valid")
    print("for fully-immersed-in-field; here ~33% iron, so p_eff is a relative metric).")
    print()
    print(f"{'win':>4} {'tid':>3} {'nMTC':>4} "
          f"{'reco R':>10} {'truth R':>10} {'R_t/R_r':>9} "
          f"{'reco p_eff':>11} {'truth p_eff':>12} "
          f"{'|dx|mean':>9} {'|dx|max':>9}"
          f" {'chi2/ndf':>9}")
    for s in summaries:
        for t in s["tracks"]:
            rfit = t.get("reco_fit")  or {}
            tfit = t.get("truth_fit") or {}
            rR = rfit.get("radius_mm", float("nan"))
            tR = tfit.get("radius_mm", float("nan"))
            ratio = tR / rR if (rR and abs(rR) > 1e-6) else float("nan")
            r_p = momentum_from_radius(rR) if rR == rR else float("nan")
            t_p = momentum_from_radius(tR) if tR == tR else float("nan")
            print(f"{s['window_id']:>4d} {t['track_id']:>3d} "
                  f"{t['n_in_mtc']:>4d} "
                  f"{rR:>10.1f} {tR:>10.1f} {ratio:>9.2f} "
                  f"{r_p:>11.3f} {t_p:>12.3f} "
                  f"{t['x_dev_mean']:>9.1f} {t['x_dev_max']:>9.1f} {t['chi2/ndf']:>9.2f}")

    if summaries:
        all_r_R = [t["reco_fit"]["radius_mm"]  for s in summaries for t in s["tracks"]
                   if t.get("reco_fit")]
        all_t_R = [t["truth_fit"]["radius_mm"] for s in summaries for t in s["tracks"]
                   if t.get("truth_fit")]
        if all_r_R and all_t_R:
            print()
            print(f"Median reco  R  = {np.median(np.abs(all_r_R)):8.1f} mm")
            print(f"Median truth R  = {np.median(np.abs(all_t_R)):8.1f} mm")
            print(f"Median |R_truth / R_reco| (over-curvature factor) = "
                  f"{np.median(np.abs(np.array(all_t_R) / np.array(all_r_R))):.2f}")
            print(f"Sign agreement (reco_curv * truth_curv > 0): "
                  f"{sum(1 for r, tt in zip(all_r_R, all_t_R) if r * tt > 0)}/{len(all_r_R)}")
    print(f"\nPlots in: {outdir}")


if __name__ == "__main__":
    main()
    # Work around a ROOT 6.38 RNTuple/RColumn teardown segfault when the
    # Python interpreter exits with cached RNTuple readers still alive.
    sys.stdout.flush(); sys.stderr.flush()
    os._exit(0)

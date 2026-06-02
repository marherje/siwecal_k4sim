#!/usr/bin/env python3
"""
SND Generic Event Display using ROOT TEve.
Detector geometry and hit display parameters are read from a JSON config file.
No hardcoded detector positions — all geometry comes from detector_config.json.

Usage:
    python event_display_eve.py \
        --hits ShipHits.root \
        --config detector_config.json \
        --window 0 \
        --source-labels "0:mixed" "1:Sig1"
"""

import argparse
import json
import math
import os
import ROOT

# Keep TColor object references alive (prevents Python GC from destroying them)
_color_refs = []

# ---------------------------------------------------------------------------
# C++ namespace for RNTuple data (avoids PyROOT GC segfault)
# ---------------------------------------------------------------------------
ROOT.gROOT.ProcessLine("""
namespace SNDEve {
    std::vector<float> gX, gY, gZ, gE;
    std::vector<int>   gSrc, gPlane;
    std::vector<float> gTsX, gTsY, gTsZ, gTsLoc0, gTsTilt;
    std::vector<int>   gTsTrackId;
}
""")


# ---------------------------------------------------------------------------
# UTILITY: nearest layer snap
# ---------------------------------------------------------------------------
def nearest_plane(hit_cm, plane_list):
    return min(plane_list, key=lambda p: abs(p - hit_cm))


# ---------------------------------------------------------------------------
# UTILITY: make TEveGeoShape box
# ---------------------------------------------------------------------------
def make_box(name, dx, dy, dz, x, y, z, color, transparency=0, stereo_deg=0.0):
    shape = ROOT.TGeoBBox(dx, dy, dz)
    gs = ROOT.TEveGeoShape(name)
    gs.SetShape(shape)
    gs.SetMainColor(color)
    gs.SetFillColor(color)
    gs.SetLineColor(color)
    gs.SetMainTransparency(transparency)
    if stereo_deg != 0.0:
        # RotateLF(1,2,phi) maps local Y → (-sin(phi), cos(phi), 0) in world.
        # stereo_deg is the physical strip angle (+5 for U, -5 for V), so we
        # negate it: phi = -stereo_deg gives local Y → (+sin(stereo_deg), cos(stereo_deg)).
        gs.RefMainTrans().RotateLF(1, 2, -stereo_deg * math.pi / 180.0)
    gs.RefMainTrans().SetPos(x, y, z)
    gs.ResetBBox()
    return gs


# ---------------------------------------------------------------------------
# UTILITY: ROOT color from RGB tuple
# ---------------------------------------------------------------------------
def rgb_to_root_color(r, g, b):
    idx = ROOT.TColor.GetFreeColorIndex()
    c = ROOT.TColor(idx, r, g, b,"")
    _color_refs.append(c)  # keep reference alive — prevents Python GC destruction
    return idx


# ---------------------------------------------------------------------------
# READ HITS FROM RNTUPLE (generic, with optional plane filter)
# ---------------------------------------------------------------------------
def read_hits(hits_file, ntuple_name, window_id, filter_dict=None, energy_field="edep"):
    """
    Read hits from an RNTuple. filter_dict applies equality filters on integer
    fields, e.g. {"plane": 0} keeps only hits where plane==0.
    energy_field selects the energy branch name: "edep" for TrackerHit3D
    measurement RNTuples, "E" for raw SimCalorimeterHit RNTuples.
    Returns xs, ys, zs, Es, srcs, planes (coords in mm, not yet converted).
    """
    filter_decls  = ""
    filter_checks = ""
    if filter_dict:
        for field, val in filter_dict.items():
            filter_decls  += f'auto vflt_{field} = reader->GetView<int>("{field}"); '
            filter_checks += f'if ((int)vflt_{field}(i) != {val}) continue; '

    ROOT.gROOT.ProcessLine(f"""
    {{
        SNDEve::gX.clear(); SNDEve::gY.clear(); SNDEve::gZ.clear();
        SNDEve::gE.clear(); SNDEve::gSrc.clear(); SNDEve::gPlane.clear();
        try {{
            auto reader = ROOT::RNTupleReader::Open("{ntuple_name}", "{hits_file}");
            auto vwin = reader->GetView<int>("window_id");
            auto vx   = reader->GetView<float>("x");
            auto vy   = reader->GetView<float>("y");
            auto vz   = reader->GetView<float>("z");
            auto vE   = reader->GetView<float>("{energy_field}");
            auto vsrc = reader->GetView<int>("source_id");
            {filter_decls}
            try {{
                auto vplane = reader->GetView<int>("plane");
                for (auto i : reader->GetEntryRange()) {{
                    if ((int)vwin(i) != {window_id}) continue;
                    {filter_checks}
                    SNDEve::gX.push_back(vx(i));
                    SNDEve::gY.push_back(vy(i));
                    SNDEve::gZ.push_back(vz(i));
                    SNDEve::gE.push_back(vE(i));
                    SNDEve::gSrc.push_back(vsrc(i));
                    SNDEve::gPlane.push_back(vplane(i));
                }}
            }} catch (...) {{
                for (auto i : reader->GetEntryRange()) {{
                    if ((int)vwin(i) != {window_id}) continue;
                    {filter_checks}
                    SNDEve::gX.push_back(vx(i));
                    SNDEve::gY.push_back(vy(i));
                    SNDEve::gZ.push_back(vz(i));
                    SNDEve::gE.push_back(vE(i));
                    SNDEve::gSrc.push_back(vsrc(i));
                    SNDEve::gPlane.push_back(-1);
                }}
            }}
        }} catch (...) {{
            printf("[reader] RNTuple '{ntuple_name}' not found or error.\\n");
        }}
        printf("[reader] %s window {window_id}: %zu hits\\n",
               "{ntuple_name}", SNDEve::gX.size());
    }}
    """)

    n      = ROOT.SNDEve.gX.size()
    xs     = [ROOT.SNDEve.gX[i]     for i in range(n)]
    ys     = [ROOT.SNDEve.gY[i]     for i in range(n)]
    zs     = [ROOT.SNDEve.gZ[i]     for i in range(n)]
    Es     = [ROOT.SNDEve.gE[i]     for i in range(n)]
    srcs   = [ROOT.SNDEve.gSrc[i]   for i in range(n)]
    planes = [ROOT.SNDEve.gPlane[i] for i in range(n)]
    return xs, ys, zs, Es, srcs, planes


# ---------------------------------------------------------------------------
# BUILD GEOMETRY from config
# ---------------------------------------------------------------------------
def build_geometry(eve, config):
    geo_list = ROOT.TEveElementList("Geometry")
    for geo in config.get("geometry", []):
        r, g, b = geo["color"]
        color   = rgb_to_root_color(r, g, b)
        transp  = geo.get("transparency", 80)
        v       = geo["voxel"]
        dx, dy, dz = v["x"], v["y"], v["z"]
        grp = ROOT.TEveElementList(geo["name"])
        for i, z in enumerate(geo["layers_z_cm"]):
            gs = make_box(f"{geo['name']}_{i}", dx, dy, dz,
                          0.0, 0.0, z, color, transparency=transp)
            grp.AddElement(gs)
        geo_list.AddElement(grp)
    eve.GetGlobalScene().AddElement(geo_list)
    print(f"[Geometry] Built {len(config.get('geometry', []))} detector groups.")


# ---------------------------------------------------------------------------
# BUILD HITS from config
# ---------------------------------------------------------------------------
def build_hits(eve, hits_file, window_id, config):
    hits_root = ROOT.TEveElementList(f"Hits window {window_id}")
    mm2cm = 0.1

    for det in config["detectors"]:
        det_name   = det["name"]
        ntuple     = det["ntuple"]
        filter_d   = det.get("filter", {})
        v          = det["voxel"]
        dx, dy, dz = v["x"], v["y"], v["z"]
        layers_z   = det["layers_z_cm"]
        z_min = det.get("z_range", {}).get("min", float("-inf"))
        z_max = det.get("z_range", {}).get("max", float("inf"))
        if len(layers_z) == 0: continue
        stereo_deg = det.get("stereo_deg", 0.0)

        xs, ys, zs, Es, srcs, planes = read_hits(
            hits_file, ntuple, window_id, filter_dict=filter_d,
            energy_field=det.get("energy_field", "edep"))

        if not xs:
            print(f"[Hits] {det_name}: no hits")
            continue

        print(f"[Hits] {det_name}: {len(xs)} hits")

        # Color all hits of this detector using the color from JSON config
        r, g, b = det["color"]
        det_color = rgb_to_root_color(r, g, b)

        det_list = ROOT.TEveElementList(det_name)
        grp = ROOT.TEveElementList(f"{det_name} hits")

        for i in range(len(xs)):
            xc = xs[i] * mm2cm
            yc = 0.0 if stereo_deg != 0.0 else ys[i] * mm2cm
            zc = zs[i] * mm2cm
            # Filter on the raw hit Z so hits from adjacent MTC stations are
            # rejected before snapping. Using snap_z here caused out-of-station
            # hits to snap to the nearest boundary layer and pass the range check.
            if not (z_min <= zc <= z_max):
                continue
            snap_z = nearest_plane(zc, layers_z)
            gs = make_box(f"{det_name}_hit_{i}",
                          dx, dy, dz,
                          xc, yc, snap_z,
                          det_color, transparency=0,
                          stereo_deg=stereo_deg)
            grp.AddElement(gs)

        det_list.AddElement(grp)
        print(f"[Hits]   {det_name}: {len(xs)} hits")
        hits_root.AddElement(det_list)

    eve.GetEventScene().AddElement(hits_root)


# ---------------------------------------------------------------------------
# GEOMETRY EXTRACTION from DD4hep (before TEveManager)
# ---------------------------------------------------------------------------
def extract_z_from_geometry(compact_xml, config):
    """
    Walk the DD4hep gGeoManager and fill layers_z_cm in-place for every
    config entry that carries a 'geo_extract' rule.
    Must be called before TEveManager.Create() — the two frameworks coexist
    because TEveGeoShape positions via RefMainTrans, independent of gGeoManager.
    Translation values from GetCurrentMatrix() are in cm (ROOT TGeo convention).
    """
    ROOT.gSystem.Load("libDDCore")
    desc = ROOT.dd4hep.Detector.getInstance()
    desc.fromXML(compact_xml)
    mgr = ROOT.gGeoManager

    rules = {}
    for entry in config.get("detectors", []) + config.get("geometry", []):
        ge = entry.get("geo_extract")
        if ge:
            rules[id(entry)] = {
                "path_contains": ge["path_contains"],
                "slice_suffix":  f"_slice_{ge['slice_index']}",
                "entry":         entry,
                "zs":            [],
            }

    def walk(node, path=""):
        vol  = node.GetVolume()
        name = vol.GetName()
        current_path = path + "/" + node.GetName()
        if "_slice_" in name:
            if mgr.cd(current_path):
                t = mgr.GetCurrentMatrix().GetTranslation()
                for r in rules.values():
                    if (r["path_contains"] in current_path
                            and name.endswith(r["slice_suffix"])):
                        r["zs"].append(round(t[2], 4))
        for i in range(node.GetNdaughters()):
            walk(node.GetDaughter(i), current_path)

    walk(mgr.GetTopNode())

    for r in rules.values():
        zs = sorted(set(r["zs"]))
        r["entry"]["layers_z_cm"] = zs
        if zs:
            # 3 cm margin covers hits anywhere inside the half-layer at the
            # station edge (~3.7 cm half-spacing). It stays within the station
            # because the gap between the last sensitive layer of one MTC station
            # and the first of the next is ~5.4 cm (1 mm gap + 50 mm Fe + 3 mm
            # inner iron), so 3 cm < 5.4 cm and no cross-station contamination
            # occurs — provided the filter is applied to the raw hit Z (zc),
            # not to the snapped Z.
            r["entry"]["z_range"] = {"min": zs[0] - 3.0, "max": zs[-1] + 3.0}
        print(f"[Geometry] {r['entry']['name']}: {len(zs)} layers extracted")


# ---------------------------------------------------------------------------
# TRACK READER (before TEveManager)
# ---------------------------------------------------------------------------
def read_track_points(hits_file, window_id):
    """Read ACTSTracks seeds from ShipHits.root before TEve is created."""
    if not os.path.exists(hits_file):
        print(f"[Tracks] ERROR: file not found: '{hits_file}'")
        raise SystemExit(1)

    ROOT.gSystem.Load("libROOTNTuple")
    ROOT.gROOT.ProcessLine(f"""
    {{
      SNDEve::gX.clear(); SNDEve::gY.clear();
      SNDEve::gZ.clear(); SNDEve::gE.clear();
      SNDEve::gSrc.clear(); SNDEve::gPlane.clear();
      try {{
        auto r       = ROOT::RNTupleReader::Open("ACTSTracks", "{hits_file}");
        auto vwin    = r->GetView<int>("window_id");
        auto vtid    = r->GetView<int>("track_id");
        auto vchi2   = r->GetView<float>("chi2");
        auto vndf    = r->GetView<int>("ndf");
        auto vseedx  = r->GetView<float>("seed_x");
        auto vseedy  = r->GetView<float>("seed_y");
        for (auto i : r->GetEntryRange()) {{
          if ((int)vwin(i) != {window_id}) continue;
          SNDEve::gX.push_back(vseedx(i));
          SNDEve::gY.push_back(vseedy(i));
          SNDEve::gZ.push_back(vchi2(i));
          SNDEve::gE.push_back((float)vndf(i));
          SNDEve::gSrc.push_back(vtid(i));
          SNDEve::gPlane.push_back(0);
        }}
      }} catch (...) {{}}
    }}
    """)

    n = ROOT.SNDEve.gX.size()
    if n == 0:
        print(f"[Tracks] No tracks in window {window_id}.")
        return []
    print(f"[Tracks] Window {window_id}: {n} track(s) found.")
    seeds = []
    for i in range(n):
        seeds.append({
            'track_id': ROOT.SNDEve.gSrc[i],
            'seed_x':   ROOT.SNDEve.gX[i],
            'seed_y':   ROOT.SNDEve.gY[i],
            'chi2':     ROOT.SNDEve.gZ[i],
            'ndf':      int(ROOT.SNDEve.gE[i]),
        })

    ROOT.gROOT.ProcessLine(f"""
    {{
      SNDEve::gZ.clear();
      try {{
        auto r   = ROOT::RNTupleReader::Open("SiTargetMeas", "{hits_file}");
        auto vw  = r->GetView<int>("window_id");
        auto vz  = r->GetView<float>("z");
        std::set<int> seen;
        for (auto i : r->GetEntryRange()) {{
          if ((int)vw(i) != {window_id}) continue;
          int iz = (int)std::round(vz(i));
          if (!seen.count(iz)) {{ seen.insert(iz); SNDEve::gZ.push_back(vz(i)); }}
        }}
      }} catch (...) {{}}
      try {{
        auto r   = ROOT::RNTupleReader::Open("SiPadMeas", "{hits_file}");
        auto vw  = r->GetView<int>("window_id");
        auto vz  = r->GetView<float>("z");
        std::set<int> seen;
        for (auto i : r->GetEntryRange()) {{
          if ((int)vw(i) != {window_id}) continue;
          int iz = (int)std::round(vz(i));
          if (!seen.count(iz)) {{ seen.insert(iz); SNDEve::gZ.push_back(vz(i)); }}
        }}
      }} catch (...) {{}}
      try {{
        auto r   = ROOT::RNTupleReader::Open("MTCSciFiMeas", "{hits_file}");
        auto vw  = r->GetView<int>("window_id");
        auto vz  = r->GetView<float>("z");
        std::set<int> seen;
        for (auto i : r->GetEntryRange()) {{
          if ((int)vw(i) != {window_id}) continue;
          int iz = (int)std::round(vz(i));
          if (!seen.count(iz)) {{ seen.insert(iz); SNDEve::gZ.push_back(vz(i)); }}
        }}
      }} catch (...) {{}}
    }}
    """)

    mm2cm = 0.1
    n_z   = ROOT.SNDEve.gZ.size()
    layer_zs_mm = sorted([ROOT.SNDEve.gZ[i] for i in range(n_z)])

    states = read_track_states(hits_file, window_id)

    track_data = []
    for seed in seeds:
        tid    = seed['track_id']
        pts_3d = states.get(tid, [])

        if len(pts_3d) >= 2:
            points = [(x, y, z, z / mm2cm) for x, y, z in pts_3d]
        elif layer_zs_mm:
            sx, sy = seed['seed_x'], seed['seed_y']
            points = [(sx*mm2cm, sy*mm2cm, z*mm2cm, z) for z in layer_zs_mm]
        else:
            continue

        if len(points) >= 2:
            track_data.append({
                'points': points,
                'chi2':   seed['chi2'],
                'ndf':    seed['ndf'],
                'seed_x': seed['seed_x'],
                'seed_y': seed['seed_y'],
            })
            print(f"[Tracks]   Track {tid}: "
                  f"chi2={seed['chi2']:.1f}  nStates={len(pts_3d)}"
                  f"  seed=({seed['seed_x']:.1f},{seed['seed_y']:.1f})mm")
    return track_data


# ---------------------------------------------------------------------------
# PER-SURFACE TRACK STATE READER
# ---------------------------------------------------------------------------
def read_track_states(hits_file, window_id):
    """Read ACTSTrackStates RNTuple; returns {track_id: [(x_cm, y_cm, z_cm), ...]}
    with stereo pairing applied for MTC SciFi U/V planes."""
    ROOT.gROOT.ProcessLine(f"""
    {{
      SNDEve::gTsX.clear(); SNDEve::gTsY.clear();
      SNDEve::gTsZ.clear(); SNDEve::gTsTrackId.clear();
      SNDEve::gTsLoc0.clear(); SNDEve::gTsTilt.clear();
      try {{
        auto r     = ROOT::RNTupleReader::Open("ACTSTrackStates", "{hits_file}");
        auto vwin  = r->GetView<int>("window_id");
        auto vtid  = r->GetView<int>("track_id");
        auto vx    = r->GetView<float>("x");
        auto vy    = r->GetView<float>("y");
        auto vz    = r->GetView<float>("z");
        auto vloc0 = r->GetView<float>("loc0");
        auto vtilt = r->GetView<float>("tilt");
        for (auto i : r->GetEntryRange()) {{
          if ((int)vwin(i) != {window_id}) continue;
          SNDEve::gTsTrackId.push_back(vtid(i));
          SNDEve::gTsX.push_back(vx(i));
          SNDEve::gTsY.push_back(vy(i));
          SNDEve::gTsZ.push_back(vz(i));
          SNDEve::gTsLoc0.push_back(vloc0(i));
          SNDEve::gTsTilt.push_back(vtilt(i));
        }}
      }} catch (...) {{}}
    }}
    """)
    n     = ROOT.SNDEve.gTsX.size()
    mm2cm = 0.1
    TAN5  = math.tan(5.0 * math.pi / 180.0)

    # Collect raw state data per track (in state_idx order, already beam-direction)
    raw = {}
    for i in range(n):
        tid = ROOT.SNDEve.gTsTrackId[i]
        raw.setdefault(tid, []).append((
            ROOT.SNDEve.gTsX[i] * mm2cm,     # x [cm] from inverse-rotation formula
            ROOT.SNDEve.gTsY[i] * mm2cm,     # y [cm] from inverse-rotation formula
            ROOT.SNDEve.gTsZ[i] * mm2cm,     # beam-Z [cm]
            float(ROOT.SNDEve.gTsLoc0[i]),   # raw ACTS eBoundLoc0 [mm]
            float(ROOT.SNDEve.gTsTilt[i]),   # stereo tilt: +sin5° U, -sin5° V, 0 non-stereo
        ))

    # Apply stereo pairing: for MTC SciFi U/V pairs in the same layer,
    # x = (loc0_U + loc0_V)/2,  y = (loc0_V - loc0_U)/(2*tan5°)
    # This avoids using loc1 (unreliable for 1D strip measurements through iron).
    by_track = {}
    for tid, states in raw.items():
        pts = []
        i = 0
        while i < len(states):
            x, y, z, loc0, tilt = states[i]
            if abs(tilt) > 0.01:
                # MTC stereo surface — look for partner in adjacent state
                if i + 1 < len(states):
                    _, _, z2, loc02, tilt2 = states[i + 1]
                    if abs(z - z2) < 0.5 and tilt * tilt2 < 0:  # same layer, opp. tilt
                        loc0_U = loc0  if tilt  > 0 else loc02   # positive tilt → V plane
                        loc0_V = loc02 if tilt  > 0 else loc0    # negative tilt → U plane
                        gx = (loc0_U + loc0_V) / 2.0 * mm2cm
                        gy = (loc0_V - loc0_U) / (2.0 * TAN5) * mm2cm
                        gz = (z + z2) / 2.0
                        pts.append((gx, gy, gz))
                        i += 2
                        continue
                # Unpaired stereo: fall back to inverse-rotation result
                pts.append((x, y, z))
            else:
                # Non-stereo (SiTarget, SiPad): inverse-rotation gives correct x,y
                pts.append((x, y, z))
            i += 1
        by_track[tid] = pts

    print(f"[Tracks] ACTSTrackStates window {window_id}: {n} raw state(s), "
          f"{sum(len(v) for v in by_track.values())} display point(s) "
          f"across {len(by_track)} track(s).")
    return by_track


# ---------------------------------------------------------------------------
# TRACK DRAWER (after TEveManager)
# ---------------------------------------------------------------------------
def draw_tracks(eve, track_data):
    if not track_data:
        return
    colors = [ROOT.kRed, ROOT.kOrange+7, ROOT.kMagenta, ROOT.kYellow+1]
    tracks_list = ROOT.TEveElementList("Reconstructed Tracks")
    for iT, td in enumerate(track_data):
        color  = colors[iT % len(colors)]
        points = td['points']
        label  = f"Track {iT}  chi2={td['chi2']:.1f}  ndf={td['ndf']}"
        line   = ROOT.TEveStraightLineSet(label)
        line.SetLineColor(color)
        line.SetLineWidth(3)
        line.SetMarkerColor(color)
        line.SetMarkerSize(0.8)
        line.SetMarkerStyle(4)
        for j in range(len(points)-1):
            x0, y0, z0, _ = points[j]
            x1, y1, z1, _ = points[j+1]
            line.AddLine(x0, y0, z0, x1, y1, z1)
        for pt in points:
            line.AddMarker(pt[0], pt[1], pt[2])
        tracks_list.AddElement(line)
    eve.GetEventScene().AddElement(tracks_list)
    print(f"[Tracks] Added {len(track_data)} track(s).")


# ---------------------------------------------------------------------------
# MAIN
# ---------------------------------------------------------------------------
def main():
    parser = argparse.ArgumentParser(description="SND Generic Event Display")
    parser.add_argument("--hits",   required=True,
                        help="Path to ShipHits.root (RNTuple file)")
    parser.add_argument("--config", default="detector_config.json",
                        help="Path to detector JSON config file")
    parser.add_argument("--geometry",
                        default=os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                             "../simulation/geometry/SND_compact.xml"),
                        help="Path to DD4hep compact XML for layer z-extraction")
    parser.add_argument("--window", type=int, default=0,
                        help="Window/event index")
    args = parser.parse_args()

    if not os.path.exists(args.config):
        print(f"[Config] ERROR: config file not found: '{args.config}'")
        raise SystemExit(1)
    with open(args.config) as f:
        config = json.load(f)
    print(f"[Config] Loaded {len(config['detectors'])} detector(s) from '{args.config}'")

    if not os.path.exists(args.geometry):
        print(f"[Geometry] ERROR: compact XML not found: '{args.geometry}'")
        raise SystemExit(1)
    extract_z_from_geometry(args.geometry, config)

    track_data = read_track_points(args.hits, args.window)

    ROOT.gErrorIgnoreLevel = ROOT.kWarning
    ROOT.gSystem.Load("libEve")
    ROOT.gSystem.Load("libEG")
    ROOT.gSystem.Load("libROOTNTuple")

    eve = ROOT.TEveManager.Create(True, "FI")
    viewer = eve.SpawnNewViewer("3D View", "")
    viewer.AddScene(eve.GetGlobalScene())
    viewer.AddScene(eve.GetEventScene())

    build_geometry(eve, config)
    build_hits(eve, args.hits, args.window, config)
    if track_data:
        draw_tracks(eve, track_data)

    eve.Redraw3D(True)
    ROOT.gSystem.ProcessEvents()

    print(f"\n[Ready] Displaying window {args.window}.")
    print("[Ready] Close the TEve window to exit.")
    ROOT.gApplication.Run(True)


if __name__ == "__main__":
    main()

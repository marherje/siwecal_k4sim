import cppyy
import ROOT
ROOT.gSystem.Load("libDDCore")

desc = ROOT.dd4hep.Detector.getInstance()
desc.fromXML("../simulation/geometry/SND_compact.xml")

mgr = ROOT.gGeoManager

z_lists = {"SiPad": [], "SiTarget_X": [], "SiTarget_Y": []}

def walk(node, path=""):
    vol = node.GetVolume()
    name = vol.GetName()
    current_path = path + "/" + node.GetName()
    shape = vol.GetShape()

    # --- Only look at slices ---
    if "_slice_" in name:
        ok = mgr.cd(current_path)
        if ok:
            t = mgr.GetCurrentMatrix().GetTranslation()

            # --- SiPad: only slice_4 ---
            if "SiPad" in current_path and "_slice_4" in name:
                print(f"{'SiPad':<12} {name:<45} x={t[2]:8.2f}")
                z_lists["SiPad"].append(round(t[2], 2))

            # --- SiTarget: slice_2 (Y) and slice_4 (Z) ---
            elif "SiTarget" in current_path:
                if "_slice_2" in name:
                    print(f"{'SiTarget_X':<12} {name:<45} x={t[2]:8.2f}")
                    z_lists["SiTarget_X"].append(round(t[2], 2))

                elif "_slice_4" in name:
                    print(f"{'SiTarget_Y':<12} {name:<45} x={t[2]:8.2f}")
                    z_lists["SiTarget_Y"].append(round(t[2], 2))

    # --- Recurse ---
    for i in range(node.GetNdaughters()):
        walk(node.GetDaughter(i), current_path)

walk(mgr.GetTopNode())

print()

var_names = {
    "SiPad": "SIPAD_Z",
    "SiTarget_X": "SITARGET_X_Z",
    "SiTarget_Y": "SITARGET_Y_Z"
}

for det, zs in z_lists.items():
    entries = [f"{z:.2f}" for z in zs]

    lines, line = [], []
    for e in entries:
        line.append(e)
        if len(", ".join(line)) > 55:
            lines.append("    " + ", ".join(line[:-1]) + ",")
            line = [line[-1]]
    if line:
        lines.append("    " + ", ".join(line) + ",")

    print(f"{var_names[det]} = [")
    print("\n".join(lines))
    print("]")
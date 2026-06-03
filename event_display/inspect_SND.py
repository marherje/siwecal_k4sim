import cppyy
import ROOT
ROOT.gSystem.Load("libDDCore")

desc = ROOT.dd4hep.Detector.getInstance()
desc.fromXML("../simulation/geometry/SND_compact.xml")

mgr = ROOT.gGeoManager

z_lists = {"SiPad": []}

def walk(node, path=""):
    vol = node.GetVolume()
    name = vol.GetName()
    current_path = path + "/" + node.GetName()
    shape = vol.GetShape()

    if "_slice_" in name:
        ok = mgr.cd(current_path)
        if ok:
            t = mgr.GetCurrentMatrix().GetTranslation()

            if "SiPad" in current_path and "_slice_4" in name:
                print(f"{'SiPad':<12} {name:<45} x={t[2]:8.2f}")
                z_lists["SiPad"].append(round(t[2], 2))

    for i in range(node.GetNdaughters()):
        walk(node.GetDaughter(i), current_path)

walk(mgr.GetTopNode())

print()

var_names = {
    "SiPad": "SIPAD_Z",
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

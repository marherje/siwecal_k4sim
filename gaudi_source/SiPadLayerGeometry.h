/**
 * SiPadLayerGeometry — where the sensitive layers are, according to DD4hep.
 *
 * One implementation of "walk the SiPad DetElement tree and find, per layer,
 * the sensitive slice and the material stack in front of it". Both ACTSGeoSvc
 * (which needs z, transverse size, thickness and material to build the ACTS
 * surfaces) and DetectorFlipper (which needs only z) read it, so there is a
 * single answer to "where is layer N" in the code.
 *
 * That mattered: the layer pitch went from 11 mm to 15 mm in July 2026, and
 * every place holding its own copy of the z table silently kept the old one.
 * The geometry is the only thing that cannot go stale, because it is what the
 * simulation actually used.
 *
 * Deliberately free of ACTS types so it can be included from algorithms that
 * have nothing to do with tracking: materials come out as (name, thickness)
 * pairs and the caller turns them into whatever it needs.
 *
 * Units: TGeo/DD4hep work in cm, everything returned here is in mm.
 */
#ifndef SIPAD_LAYER_GEOMETRY_H
#define SIPAD_LAYER_GEOMETRY_H

#include "DD4hep/DetElement.h"
#include "DD4hep/Detector.h"
#include "DD4hep/Volumes.h"

#include "TGeoBBox.h"
#include "TGeoMatrix.h"
#include "TGeoNode.h"

#include <functional>
#include <map>
#include <stdexcept>
#include <string>
#include <vector>

namespace sipad {

/// One non-air slice of a layer: what it is made of and how thick it is.
struct MaterialSlice {
  std::string material;
  double      thicknessMm = 0.0;
};

/// The sensitive plane of one layer, plus the material stack of that layer.
struct LayerPlane {
  int    layer         = -1;   ///< DetElement id = layer number
  double z             = 0.0;  ///< world z of the sensitive slice centre [mm]
  double halfX         = 0.0;  ///< transverse half-size [mm]
  double halfY         = 0.0;  ///< transverse half-size [mm]
  double halfThickness = 0.0;  ///< half-thickness along z [mm]
  std::vector<MaterialSlice> slices;  ///< every non-air slice, in slice order
};

namespace detail {

/// True when the volume or any descendant carries the sensitive-detector flag.
/// With the tiled segmentation the flag sits on the <slice>_wafer_pads volume
/// (detector_plugin/SiPadDetector.cpp, buildWafers), not on the slice
/// container, which is plain gap material — so asking the container alone
/// finds nothing.
inline bool subtreeIsSensitive(TGeoVolume* v) {
  if (!v) return false;
  if (dd4hep::Volume(v).isSensitive()) return true;
  for (int i = 0; i < v->GetNdaughters(); ++i)
    if (subtreeIsSensitive(v->GetNode(i)->GetVolume())) return true;
  return false;
}

inline TGeoVolume* firstSensitiveVolume(TGeoVolume* v) {
  if (!v) return nullptr;
  if (dd4hep::Volume(v).isSensitive()) return v;
  for (int i = 0; i < v->GetNdaughters(); ++i)
    if (TGeoVolume* s = firstSensitiveVolume(v->GetNode(i)->GetVolume()))
      return s;
  return nullptr;
}

}  // namespace detail

/**
 * One LayerPlane per layer of @p detName, ordered by layer number.
 *
 * The tree is walked by DetElement rather than by matching TGeo volume names:
 * the SiPadDetector plugin registers one DetElement per layer (id = layer
 * number) and one child per slice (id = slice number), so the sensitive slice
 * is found by asking DD4hep about sensitivity. Name matching is what broke
 * before, when a segmentation change moved the silicon from slice 4 to slice 5
 * and the hardcoded "_slice_4" match quietly started returning air.
 *
 * @throws std::runtime_error if the detector is absent, or if any layer does
 *         not have exactly one sensitive slice — a geometry the surface model
 *         cannot represent, which is worth stopping for rather than guessing.
 */
inline std::vector<LayerPlane> sensitiveLayers(dd4hep::Detector& desc,
                                               const std::string& detName = "SiPad") {
  dd4hep::DetElement sdet = desc.detector(detName);
  if (!sdet.isValid())
    throw std::runtime_error("SiPadLayerGeometry: detector '" + detName +
                             "' not found in the compact geometry");

  // children() is keyed by name; re-key by DetElement id (= layer number) so
  // the result is ordered by layer and not by however the names happen to sort.
  std::map<int, dd4hep::DetElement> layersById;
  for (const auto& [name, layerDe] : sdet.children())
    layersById[layerDe.id()] = layerDe;

  std::vector<LayerPlane> planes;
  planes.reserve(layersById.size());

  for (const auto& [layerId, layerDe] : layersById) {
    std::map<int, dd4hep::DetElement> slicesById;
    for (const auto& [name, sliceDe] : layerDe.children())
      slicesById[sliceDe.id()] = sliceDe;

    LayerPlane plane;
    plane.layer   = layerId;
    int nSensitive = 0;

    for (const auto& [sliceId, sliceDe] : slicesById) {
      TGeoVolume* vol = sliceDe.volume().ptr();
      TGeoBBox*   box = vol ? dynamic_cast<TGeoBBox*>(vol->GetShape()) : nullptr;
      if (!box) continue;

      const bool sensitive = detail::subtreeIsSensitive(vol);
      if (sensitive) {
        ++nSensitive;
        const TGeoHMatrix& w = sliceDe.nominal().worldTransformation();
        plane.z             = w.GetTranslation()[2] * 10.0;
        plane.halfX         = box->GetDX() * 10.0;
        plane.halfY         = box->GetDY() * 10.0;
        plane.halfThickness = box->GetDZ() * 10.0;
      }

      // Material budget: every non-air slice contributes. For a tiled
      // sensitive slice the container is gap material (air), so take the
      // material of the sensitive volume itself — the silicon wafer.
      std::string matName = dd4hep::Volume(vol).material().name();
      if (sensitive) {
        if (TGeoVolume* sv = detail::firstSensitiveVolume(vol))
          matName = dd4hep::Volume(sv).material().name();
      }
      if (matName != "Air" && matName != "Vacuum")
        plane.slices.push_back({matName, 2.0 * box->GetDZ() * 10.0});
    }

    if (nSensitive != 1)
      throw std::runtime_error(
          "SiPadLayerGeometry: layer " + std::to_string(layerId) + " of " +
          detName + " has " + std::to_string(nSensitive) +
          " sensitive slices (expected 1) — geometry and surface model out of sync");

    planes.push_back(std::move(plane));
  }
  return planes;
}

/// Just the per-layer z [mm], ordered by layer number.
inline std::vector<double> layerZ(dd4hep::Detector& desc,
                                  const std::string& detName = "SiPad") {
  std::vector<double> z;
  for (const auto& p : sensitiveLayers(desc, detName)) z.push_back(p.z);
  return z;
}

}  // namespace sipad

#endif  // SIPAD_LAYER_GEOMETRY_H

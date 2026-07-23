//====================================================================
//  DD4hep detector plugin for SiPadDetector
//  Envelope placed via createPlacedEnvelope (DD4hep framework requirement).
//  Envelope z extent computed from layer thicknesses; dim.z() is not used
//  for layer positioning.
//
//  Sensitive slices are tiled: a plane of inert gap material holds an array of
//  silicon sensors, each of them a dead rim around a sensitive pad array. The
//  dead regions of a layer -- the rim, the clearance between sensors and the
//  outer border -- are physical. See buildWafers() below.
//====================================================================
#include "DD4hep/DetFactoryHelper.h"
#include "XML/Utilities.h"

#include <cmath>
#include <stdexcept>
#include <tuple>

using namespace dd4hep;
using namespace dd4hep::detail;

namespace {

  /// Wafer layout of a sensitive slice, read from the <tiling/> element.
  ///
  /// Two levels: a plane holds asus_x * asus_y ASUs on a fixed pitch, and each
  /// ASU holds n_x * n_y wafers separated by the guard-ring gap. The dead region
  /// between two ASUs is not a parameter: it follows from the ASU pitch minus
  /// the span the wafers actually occupy.
  struct Tiling {
    int         asus_x  = 1,   asus_y  = 1;   // ASUs per plane
    double      pitch_x = 0.0, pitch_y = 0.0; // ASU pitch (ignored for a single ASU)
    int         n_x     = 0,   n_y     = 0;   // wafers per ASU
    double      wafer_x = 0.0, wafer_y = 0.0; // physical sensor size
    double      marg_x  = 0.0, marg_y  = 0.0; // inactive rim of a sensor
    double      gap_x   = 0.0, gap_y   = 0.0; // extra clearance between sensors
    std::string material = "Air";             // fills the gaps and the border

    /// Sensitive area of one sensor: the sensor minus its rim on both sides.
    double activeX() const { return wafer_x - 2.0 * marg_x; }
    double activeY() const { return wafer_y - 2.0 * marg_y; }
    /// Span of the wafers of one ASU.
    double asuSpanX() const { return n_x * wafer_x + (n_x - 1) * gap_x; }
    double asuSpanY() const { return n_y * wafer_y + (n_y - 1) * gap_y; }
    /// Span of the whole wafer array, edge to edge.
    double spanX() const { return (asus_x - 1) * pitch_x + asuSpanX(); }
    double spanY() const { return (asus_y - 1) * pitch_y + asuSpanY(); }
    bool   valid() const {
      return n_x > 0 && n_y > 0 && activeX() > 0.0 && activeY() > 0.0;
    }
  };

  /// Populate a (non-sensitive) plane with the array of silicon sensors.
  ///
  /// A sensor is silicon throughout, but only its central 'activeX x activeY'
  /// is a diode: the rim of 'marg_x' / 'marg_y' around it is dead. It is
  /// therefore built as two nested boxes, a non-sensitive sensor body holding a
  /// sensitive pad array, so that the rim keeps the correct material budget
  /// while producing no hit.
  ///
  /// The other dead regions need no volume at all: they are simply the part of
  /// the plane that no sensor covers, i.e. the 'gap_material' the plane is made
  /// of -- the optional clearance between sensors, and the outer border.
  ///
  /// The 'wafer' physVolID sits on the sensor body, so it is picked up by every
  /// hit in the pad array inside it. The readout segmentation consequently works
  /// in the pad-array frame, where the pads form a plain uniform grid on the
  /// pixel pitch, and uses that same physVolID to pick the offset giving global
  /// pad indices.
  ///
  /// The physVolID numbers the wafers over the whole plane, row by row, ignoring
  /// the ASU they belong to: id = gy * (asus_x * n_x) + gx. That way the readout
  /// only needs the wafer's global column/row to know its offset, and a single
  /// ASU reproduces the plain 0..n_x*n_y-1 numbering.
  void buildWafers(Detector& description, Volume& plane, const std::string& name,
                   Material silicon, double thickness, SensitiveDetector& sens,
                   const std::string& vis, const Tiling& t) {
    // The sensitive pad array. When there is no rim it is the whole sensor and
    // no body volume is built, which keeps the plain single-box geometry.
    Volume pads(name + "_wafer_pads",
                Box(t.activeX() / 2.0, t.activeY() / 2.0, thickness / 2.0), silicon);
    pads.setSensitiveDetector(sens);
    pads.setVisAttributes(description, vis);

    const bool rim = t.marg_x > 0.0 || t.marg_y > 0.0;
    Volume sensor = pads;
    if (rim) {
      sensor = Volume(name + "_wafer",
                      Box(t.wafer_x / 2.0, t.wafer_y / 2.0, thickness / 2.0), silicon);
      sensor.setVisAttributes(description, vis);
      sensor.placeVolume(pads, Position(0.0, 0.0, 0.0));
    }

    const int cols = t.asus_x * t.n_x;

    for (int ay = 0; ay < t.asus_y; ++ay) {
      for (int ax = 0; ax < t.asus_x; ++ax) {
        double asu_x = (ax - (t.asus_x - 1) / 2.0) * t.pitch_x;
        double asu_y = (ay - (t.asus_y - 1) / 2.0) * t.pitch_y;

        for (int wy = 0; wy < t.n_y; ++wy) {
          for (int wx = 0; wx < t.n_x; ++wx) {
            double x = asu_x + (wx - (t.n_x - 1) / 2.0) * (t.wafer_x + t.gap_x);
            double y = asu_y + (wy - (t.n_y - 1) / 2.0) * (t.wafer_y + t.gap_y);
            PlacedVolume pv = plane.placeVolume(sensor, Position(x, y, 0.0));
            pv.addPhysVolID("wafer", (ay * t.n_y + wy) * cols + (ax * t.n_x + wx));
          }
        }
      }
    }
  }

} // namespace

dd4hep::Ref_t create_SiPadDetector(dd4hep::Detector& description,
                             xml_h element,
                             dd4hep::SensitiveDetector sens) {

  xml_det_t   x_det    = element;
  std::string det_name = x_det.nameStr();
  int         det_id   = x_det.id();
  DetElement  sdet(det_name, det_id);

  // DD4hep framework requirement: createPlacedEnvelope reads <envelope> from XML.
  Volume envelope = xml::createPlacedEnvelope(description, element, sdet);
  xml::setDetectorTypeFlag(element, sdet);

  if (description.buildType() == BUILD_ENVELOPE) return sdet;

  Material air = description.air();
  sens.setType("calorimeter");

  xml_dim_t dim  = x_det.dimensions();
  double    env_w = dim.x();
  double    env_h = dim.y();

  // -- Transverse tiling of the sensitive slices -----------------------------
  Tiling tiling;
  xml_h x_tiling = x_det.child(_Unicode(tiling), false);
  if (x_tiling.ptr()) {
    xml_comp_t t = x_tiling;
    if (t.hasAttr(_Unicode(asus_x))) tiling.asus_x = t.attr<int>(_Unicode(asus_x));
    if (t.hasAttr(_Unicode(asus_y))) tiling.asus_y = t.attr<int>(_Unicode(asus_y));
    if (tiling.asus_x > 1) tiling.pitch_x = t.attr<double>(_Unicode(asu_pitch_x));
    if (tiling.asus_y > 1) tiling.pitch_y = t.attr<double>(_Unicode(asu_pitch_y));
    tiling.n_x      = t.attr<int>(_Unicode(wafers_x));
    tiling.n_y      = t.attr<int>(_Unicode(wafers_y));
    tiling.wafer_x  = t.attr<double>(_Unicode(wafer_size_x));
    tiling.wafer_y  = t.attr<double>(_Unicode(wafer_size_y));
    // Optional: a sensor with no inactive rim is entirely sensitive.
    if (t.hasAttr(_Unicode(margin_x))) tiling.marg_x = t.attr<double>(_Unicode(margin_x));
    if (t.hasAttr(_Unicode(margin_y))) tiling.marg_y = t.attr<double>(_Unicode(margin_y));
    tiling.gap_x    = t.attr<double>(_Unicode(gap_x));
    tiling.gap_y    = t.attr<double>(_Unicode(gap_y));
    tiling.material = t.attr<std::string>(_Unicode(gap_material));

    for (auto [margin, size, axis] : {std::tuple{tiling.marg_x, tiling.wafer_x, 'x'},
                                      std::tuple{tiling.marg_y, tiling.wafer_y, 'y'}})
      if (margin < 0.0 || 2.0 * margin >= size)
        throw std::runtime_error(det_name + ": sensor margin along " + axis +
                                 " leaves no active area");

    // The wafer array must fit inside the plane; whatever is left over becomes
    // the inactive border. ASUs must not overlap either.
    for (auto [span, extent, axis] : {std::tuple{tiling.spanX(), env_w, 'x'},
                                      std::tuple{tiling.spanY(), env_h, 'y'}})
      if (span > extent + 1e-9)
        throw std::runtime_error(det_name + ": wafer array along " + axis +
                                 " is wider than the detector plane");

    for (auto [n, pitch, span, axis] : {std::tuple{tiling.asus_x, tiling.pitch_x, tiling.asuSpanX(), 'x'},
                                        std::tuple{tiling.asus_y, tiling.pitch_y, tiling.asuSpanY(), 'y'}})
      if (n > 1 && pitch < span - 1e-9)
        throw std::runtime_error(det_name + ": ASU pitch along " + axis +
                                 " is smaller than the wafers it must hold");
  }

  // -- First pass: compute total Z from slice thicknesses --------------------
  double total_z = 0.0;
  for (xml_coll_t lColl(x_det, _U(layer)); lColl; ++lColl) {
    xml_comp_t x_layer = lColl;
    int repeat = x_layer.hasAttr(_U(repeat)) ? x_layer.attr<int>(_U(repeat)) : 1;
    double layer_thick = 0.0;
    for (xml_coll_t sColl(x_layer, _U(slice)); sColl; ++sColl)
      layer_thick += xml_comp_t(sColl).attr<double>(_U(thickness));
    total_z += repeat * layer_thick;
  }

  // -- Place layers along Z (beam axis) using computed total_z ---------------
  int    layer_num = 0;
  double cur_z     = -total_z / 2.0;

  for (xml_coll_t lColl(x_det, _U(layer)); lColl; ++lColl) {
    xml_comp_t x_layer = lColl;
    int repeat = x_layer.hasAttr(_U(repeat)) ? x_layer.attr<int>(_U(repeat)) : 1;

    double layer_thick = 0.0;
    for (xml_coll_t sColl(x_layer, _U(slice)); sColl; ++sColl)
      layer_thick += xml_comp_t(sColl).attr<double>(_U(thickness));

    for (int j = 0; j < repeat; j++, ++layer_num) {
      std::string  layer_name = det_name + "_layer_" + std::to_string(layer_num);
      // Z=beam: layer box wide in X and Y (transverse), thin along Z (beam).
      Volume       layer_vol(layer_name, Box(env_w / 2.0, env_h / 2.0, layer_thick / 2.0), air);
      layer_vol.setVisAttributes(description, x_layer.visStr());

      double       layer_center_z = cur_z + layer_thick / 2.0;
      PlacedVolume layer_pv = envelope.placeVolume(layer_vol, Position(0.0, 0.0, layer_center_z));
      layer_pv.addPhysVolID("layer", layer_num);

      DetElement layer_de(sdet, layer_name, layer_num);
      layer_de.setPlacement(layer_pv);

      // -- Slices ------------------------------------------------------------
      double local_z   = -layer_thick / 2.0;
      int    slice_num = 0;

      for (xml_coll_t sColl(x_layer, _U(slice)); sColl; ++sColl, ++slice_num) {
        xml_comp_t  x_slice     = sColl;
        double      slice_thick = x_slice.thickness();
        Material    slice_mat   = description.material(x_slice.materialStr());
        std::string slice_name  = layer_name + "_slice_" + std::to_string(slice_num);

        // A tiled sensitive slice is a plane of gap material holding the wafers;
        // any other slice is a plain block of its own material.
        const bool tiled     = x_slice.isSensitive() && tiling.valid();
        Material   plane_mat = tiled ? description.material(tiling.material) : slice_mat;

        // Z=beam: slice box wide in X and Y (transverse), thin along Z (beam).
        Volume slice_vol(slice_name, Box(env_w / 2.0, env_h / 2.0, slice_thick / 2.0), plane_mat);
        slice_vol.setAttributes(description, x_slice.regionStr(), x_slice.limitsStr(),
                                tiled ? "Invisible" : x_slice.visStr());

        if (tiled)
          buildWafers(description, slice_vol, slice_name, slice_mat, slice_thick, sens,
                      x_slice.visStr(), tiling);
        else if (x_slice.isSensitive())
          slice_vol.setSensitiveDetector(sens);

        double       sl_center_z = local_z + slice_thick / 2.0;
        PlacedVolume slice_pv = layer_vol.placeVolume(slice_vol, Position(0.0, 0.0, sl_center_z));
        slice_pv.addPhysVolID("slice", slice_num);

        DetElement slice_de(layer_de, slice_name, slice_num);
        slice_de.setPlacement(slice_pv);

        local_z += slice_thick;
      }

      cur_z += layer_thick;
    }
  }

  return sdet;
}

DECLARE_DETELEMENT(SiPadDetector, create_SiPadDetector)

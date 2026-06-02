//====================================================================
//  DD4hep detector plugin for SiPadDetector
//  Envelope placed via createPlacedEnvelope (DD4hep framework requirement).
//  Envelope z extent computed from layer thicknesses; dim.z() is not used
//  for layer positioning.
//====================================================================
#include "DD4hep/DetFactoryHelper.h"
#include "XML/Utilities.h"

using namespace dd4hep;
using namespace dd4hep::detail;

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

        // Z=beam: slice box wide in X and Y (transverse), thin along Z (beam).
        Volume slice_vol(slice_name, Box(env_w / 2.0, env_h / 2.0, slice_thick / 2.0), slice_mat);
        slice_vol.setAttributes(description, x_slice.regionStr(), x_slice.limitsStr(), x_slice.visStr());

        if (x_slice.isSensitive())
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

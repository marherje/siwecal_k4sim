#pragma once

#include "ISNDGeoSvc.h"
#include "GaudiKernel/Service.h"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Geometry/TrackingGeometry.hpp"
#include "Acts/MagneticField/MagneticFieldContext.hpp"
#include "Acts/Geometry/DetectorElementBase.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include <map>
#include <memory>
#include <string>
#include <tuple>
#include <unordered_map>
#include <vector>

class ACTSGeoSvc : public extends<Service, ISNDGeoSvc> {
public:
  ACTSGeoSvc(const std::string& name, ISvcLocator* svcLoc);

  StatusCode initialize() override;
  StatusCode finalize()   override;

  // IActsGeoSvc interface
  const Acts::TrackingGeometry& trackingGeometry() const override;

  // ISNDGeoSvc extensions
  const std::vector<const Acts::Surface*>& allSurfaces() const override;
  const Acts::GeometryContext& geometryContext() const override;
  const Acts::Surface* surfaceByAddress(
      int detID, int station, int layer, int plane) const override;

private:
  Gaudi::Property<std::string> m_compactFile{
      this, "CompactFile", "",
      "Path to DD4hep compact XML"};

  std::shared_ptr<const Acts::TrackingGeometry> m_trackingGeometry;
  std::vector<const Acts::Surface*>             m_allSurfaces;
  Acts::GeometryContext                         m_gctx;
  // Detector elements must outlive all surfaces (surfaces hold raw pointers back).
  std::vector<std::shared_ptr<Acts::DetectorElementBase>> m_detectorElements;

  // key = (detID, station, layerInDet, plane)
  std::map<std::tuple<int,int,int,int>, const Acts::Surface*> m_surfaceByAddressMap;
};

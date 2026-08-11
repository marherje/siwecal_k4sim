cmake_minimum_required(VERSION 3.14)
project(SND_sim)

find_package(DD4hep REQUIRED COMPONENTS DDCore DDG4)
find_package(Gaudi REQUIRED)
find_package(k4FWCore REQUIRED)
find_package(EDM4HEP REQUIRED)
find_package(ROOT REQUIRED COMPONENTS ROOTNTuple Hist)
find_package(Acts 44.3.0 REQUIRED COMPONENTS Core PluginDD4hep)
find_package(k4ActsTracking REQUIRED)

dd4hep_add_plugin(SiPadDetector_plugin SOURCES detector_plugin/SiPadDetector.cpp)
target_link_libraries(SiPadDetector_plugin DD4hep::DDCore)

gaudi_add_module(SND_reco
  SOURCES
    gaudi_source/GeV2MIPConversion.cpp
    gaudi_source/BasicDigitizer.cpp
    gaudi_source/ContributionInspector.cpp
    gaudi_source/EventWindowSplitter.cpp
    gaudi_source/EventShuffler.cpp
    gaudi_source/EDM4HEP2RNTuple.cpp
    gaudi_source/ACTSGeoSvc.cpp
    gaudi_source/ACTSProtoTracker.cpp
    gaudi_source/SiPadMeasConverter.cpp
    gaudi_source/MIPExtractor.cpp
    gaudi_source/RealDigitizer.cpp
    gaudi_source/DetectorFlipper.cpp
    gaudi_source/ChannelMapper.cpp
    gaudi_source/ShowerTagger.cpp
  LINK
    Gaudi::GaudiKernel
    k4FWCore::k4FWCore
    EDM4HEP::edm4hep
    podio::podioRootIO
    DD4hep::DDCore
    DD4hep::DDRec
    ROOT::ROOTNTuple
    ROOT::Hist
    Acts::Core
    Acts::PluginDD4hep
    k4ActsTracking::k4ActsTracking
)

 set_target_properties(SND_reco PROPERTIES
    INSTALL_RPATH_USE_LINK_PATH TRUE
    BUILD_RPATH_USE_ORIGIN      TRUE
)
# ACTS headers. Derived from the Acts::Core target that find_package(Acts) sets
# up, not hardcoded to one key4hep release's install path -- the release (and so
# the hashed acts/<version>-<hash>/include directory) is pinned in
# .key4hep-release and changes between stacks.
get_target_property(_acts_include_dirs Acts::Core INTERFACE_INCLUDE_DIRECTORIES)
if(_acts_include_dirs)
    target_include_directories(SND_reco PRIVATE ${_acts_include_dirs})
endif()

install(TARGETS SiPadDetector_plugin LIBRARY DESTINATION lib)
install(TARGETS SND_reco LIBRARY DESTINATION lib)

cmake_minimum_required(VERSION 3.16)
project(SiWECALTimeDigiExample LANGUAGES CXX)

find_package(EDM4HEP REQUIRED)
find_package(podio REQUIRED)
find_package(ROOT REQUIRED COMPONENTS Core Graf Hist Gpad)

add_executable(plot_sipad_shaping plot_sipad_shaping.cpp)
target_compile_features(plot_sipad_shaping PRIVATE cxx_std_17)
target_include_directories(plot_sipad_shaping PRIVATE
  ${CMAKE_CURRENT_SOURCE_DIR}/../../include
)
target_link_libraries(plot_sipad_shaping PRIVATE
  EDM4HEP::edm4hep
  podio::podioRootIO
  ROOT::Core
  ROOT::Graf
  ROOT::Hist
  ROOT::Gpad
)

add_executable(doublehit doublehit.cpp)
target_compile_features(doublehit PRIVATE cxx_std_17)
target_include_directories(doublehit PRIVATE
  ${CMAKE_CURRENT_SOURCE_DIR}/../../include
)
target_link_libraries(doublehit PRIVATE
  EDM4HEP::edm4hep
  podio::podioRootIO
  ROOT::Core
  ROOT::Graf
  ROOT::Hist
  ROOT::Gpad
)

add_executable(compare_sipad_digitized compare_sipad_digitized.cpp)
target_compile_features(compare_sipad_digitized PRIVATE cxx_std_17)
target_include_directories(compare_sipad_digitized PRIVATE
  ${CMAKE_CURRENT_SOURCE_DIR}/../../include
)
target_link_libraries(compare_sipad_digitized PRIVATE
  EDM4HEP::edm4hep
  podio::podioRootIO
)

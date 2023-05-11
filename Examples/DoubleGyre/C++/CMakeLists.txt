cmake_minimum_required(VERSION 3.13)
project (DoubleGyre C CXX)

include (GNUInstallDirs)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/${CMAKE_INSTALL_BINDIR}")
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/${CMAKE_INSTALL_LIBDIR}")

# since we use C++11 in this example.
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

set(_SRCS
  main.cxx
  double_gyre.cxx
)

find_package(MPI COMPONENTS C CXX)
if (NOT MPI_FOUND)
  message(STATUS
    "Skipping example: ${CMAKE_PROJECT_NAME} requires MPI.")
  return ()
endif ()

set(INSITU None CACHE STRING "Enable in-situ support")
set_property(CACHE INSITU PROPERTY STRINGS None Catalyst Ascent)

if(INSITU STREQUAL "None")
  add_executable(double_gyre ${_SRCS})

elseif(INSITU STREQUAL "Catalyst")
  add_executable(double_gyre_catalyst ${_SRCS})
  find_package(catalyst REQUIRED PATHS "${EBROOTCATALYST}/lib/cmake/catalyst-2.0")
  target_compile_definitions(double_gyre_catalyst PUBLIC USE_CATALYST=1)
  target_link_libraries(double_gyre_catalyst PRIVATE catalyst::catalyst)
elseif(INSITU STREQUAL "Ascent")
  find_package(Ascent REQUIRED PATHS "/apps/daint/UES/Ascent/ascent-install/lib/cmake/ascent")
  add_executable(double_gyre_ascent ${_SRCS} double_gyre_ascent.cxx)
  target_compile_definitions(double_gyre_ascent PRIVATE USE_ASCENT=1)
  target_link_libraries(double_gyre_ascent PRIVATE ascent::ascent )
endif()

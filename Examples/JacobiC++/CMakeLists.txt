cmake_minimum_required(VERSION 3.13)
project (PJACOBIV2 CXX)

include (GNUInstallDirs)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/${CMAKE_INSTALL_BINDIR}")
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/${CMAKE_INSTALL_LIBDIR}")

# since we use C++11 in this example.
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Since this example uses MPI, find and link against it.
find_package(MPI COMPONENTS CXX)
if (NOT MPI_FOUND)
  message(STATUS
    "Skipping example: ${CMAKE_PROJECT_NAME} requires MPI.")
  return ()
endif ()

set(C_SRCS
  solvers.cxx
  JacobiDriver.cxx
)

add_executable(pjacobi ${C_SRCS})
target_link_libraries(pjacobi PRIVATE MPI::MPI_CXX m)

set(INSITU None CACHE STRING "Enable in-situ support")
set_property(CACHE INSITU PROPERTY STRINGS None Catalyst Ascent)

if(INSITU STREQUAL "Catalyst")
  find_package(catalyst REQUIRED PATHS)
  target_compile_definitions(pjacobi
    PRIVATE
      USE_CATALYST=1)
  target_link_libraries(pjacobi PRIVATE catalyst::catalyst)
elseif(INSITU STREQUAL "Ascent")
  find_package(Ascent REQUIRED PATHS "/apps/daint/UES/Ascent/ascent-install/lib/cmake/ascent")
  target_compile_definitions(pjacobi
    PRIVATE
      USE_ASCENT=1)
  target_link_libraries(pjacobi PRIVATE ascent::ascent_mpi)
endif()

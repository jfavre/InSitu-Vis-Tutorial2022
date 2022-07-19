cmake_minimum_required(VERSION 3.13)
project (PJACOBIV2 C CXX)

include (GNUInstallDirs)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/${CMAKE_INSTALL_BINDIR}")
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/${CMAKE_INSTALL_LIBDIR}")

# since we use C++11 in this example.
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Since this example uses MPI, find and link against it.
find_package(MPI COMPONENTS C CXX)
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
target_link_libraries(pjacobi PRIVATE MPI::MPI_C MPI::MPI_CXX m)

option(USE_CATALYST "Use ParaView's Catalyst" ON)

option(INSITU_COUPLING "create an executable linked with in-situ library" ON)
if(INSITU_COUPLING)
  if(USE_CATALYST)
    find_package(catalyst REQUIRED PATHS "${ParaView_DIR}/catalyst")
    # we only need C++ compilers if we're building the Adaptor
    enable_language(CXX)

  add_executable(pjacobi_catalystV2 ${C_SRCS})
  target_compile_definitions(pjacobi_catalystV2
    PRIVATE
      "PARAVIEW_IMPL_DIR=\"${ParaView_CATALYST_DIR}\""
      USE_CATALYST=1)
  target_link_libraries(pjacobi_catalystV2 PRIVATE MPI::MPI_C MPI::MPI_CXX m catalyst::catalyst)

  endif(USE_CATALYST)
endif(INSITU_COUPLING)

project(Heat-diffusion-adios)
cmake_minimum_required(VERSION 3.26)

find_package(Fides REQUIRED)
find_package(ADIOS2 REQUIRED)

if (ADIOS2_VERSION VERSION_GREATER_EQUAL "2.8.0")
  set(ADIOS_DEFS FIDES_ADIOS_HAS_RANDOM_ACCESS)
endif()

# need to continue to support v1.9 for awhile
find_package(VTKm REQUIRED)

if (VTKm_ENABLE_CUDA)
  enable_language(CUDA)
endif()

if(ADIOS2_HAVE_MPI)
  set(FIDES_USE_MPI TRUE)
  set(MPI_DEFS FIDES_USE_MPI)
  set(MPI_LIBS MPI::MPI_CXX)
endif()

if (FIDES_USE_MPI)
  add_executable(fides-sst-reader fides-sst-reader.cpp)
  vtkm_add_target_information(fides-sst-reader
                            DROP_UNUSED_SYMBOLS
                            MODIFY_CUDA_FLAGS
                            DEVICE_SOURCES fides-sst-reader.cpp)
  if (VTKm_VERSION VERSION_EQUAL "1.9")
    set(VTKM_LIBS vtkm_cont vtkm_rendering vtkm_filter)
  elseif(VTKm_VERSION VERSION_GREATER_EQUAL "2.0")
    set(VTKM_LIBS vtkm::cont vtkm::rendering vtkm::filter)
  endif()
  target_link_libraries(fides-sst-reader PUBLIC
    adios2::adios2 ${VTKM_LIBS} MPI::MPI_C MPI::MPI_CXX fides)
endif()


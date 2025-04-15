cmake_minimum_required(VERSION 3.12)
cmake_policy(VERSION 3.12)

enable_language(CXX)

# Set project name
project(MPulse VERSION 1.0.0)

# Set C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Set build type if not specified
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

# Set compiler flags
set(CMAKE_CXX_FLAGS "-Wall")
set(CMAKE_CXX_FLAGS_DEBUG "-g -O0")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")

# MPI setup
find_package(MPI REQUIRED)
set(CMAKE_CXX_COMPILER ${MPI_CXX_COMPILER})

# Define dimensions to build for
set(DIMENSIONS 1 2 3)

# Define include directories
set(INCLUDE_DIRS
    /usr/local/include
    /home/terence411/resources/include
    /home/terence411/spack/opt/spack/linux-ubuntu22.04-skylake/gcc-12.3.0/boost-1.82.0-3zvrwkhbsxoaivfmmy2gonv4qwdn36fb/include
    /home/terence411/spack/opt/spack/linux-ubuntu22.04-skylake/gcc-12.3.0/hdf5-1.14.5-6ibf3iy452splddw7bjwtnmu3ayj6q5t/include
    /home/terence411/spack/opt/spack/linux-ubuntu22.04-skylake/gcc-12.3.0/kokkos-4.4.01-pytybmwfndyyu6x2f5uqo5ucjgfojxzp/include
)

# Define libraries
set(LIB_DIRS
    /usr/local/lib
    /home/terence411/resources/lib
    /home/terence411/spack/opt/spack/linux-ubuntu22.04-skylake/gcc-12.3.0/hdf5-1.14.5-6ibf3iy452splddw7bjwtnmu3ayj6q5t/lib
    /home/terence411/spack/opt/spack/linux-ubuntu22.04-skylake/gcc-12.3.0/kokkos-4.4.01-pytybmwfndyyu6x2f5uqo5ucjgfojxzp/lib
)

# Set library paths for runtime
set(CMAKE_INSTALL_RPATH ${LIB_DIRS})
set(CMAKE_BUILD_WITH_INSTALL_RPATH TRUE)

# Find required packages
find_library(HDF5_LIBRARY hdf5 PATHS ${LIB_DIRS})
find_library(SCHNEK_LIBRARY schnek PATHS ${LIB_DIRS})
find_library(FFTW_LIBRARY fftw3 PATHS ${LIB_DIRS})
find_library(KOKKOS_CORE_LIBRARY kokkoscore PATHS ${LIB_DIRS})
find_library(KOKKOS_CONTAINERS_LIBRARY kokkoscontainers PATHS ${LIB_DIRS})
find_library(KOKKOS_SIMD_LIBRARY kokkossimd PATHS ${LIB_DIRS})

# Define source files
set(SOURCE_FILES
    ${CMAKE_CURRENT_SOURCE_DIR}/src/*.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/huerto/electromagnetics/current.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/huerto/electromagnetics/em_fields.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/huerto/electromagnetics/fdtd/fdtd_plain.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/huerto/electromagnetics/source/border.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/huerto/electromagnetics/source/incsource.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/huerto/electromagnetics/source/plane_wave.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/huerto/electromagnetics/source/beam.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/huerto/electromagnetics/pml/cpml_border.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/huerto/maths/random.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/huerto/maths/functions/core.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/huerto/simulation/task.cpp
)

# Use glob to expand src/*.cpp
file(GLOB SRC_CPP_FILES ${CMAKE_CURRENT_SOURCE_DIR}/src/*.cpp)

# Create a list of all source files (excluding the glob pattern)
set(ALL_SOURCE_FILES
    ${SRC_CPP_FILES}
    ${CMAKE_CURRENT_SOURCE_DIR}/huerto/electromagnetics/current.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/huerto/electromagnetics/em_fields.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/huerto/electromagnetics/fdtd/fdtd_plain.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/huerto/electromagnetics/source/border.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/huerto/electromagnetics/source/incsource.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/huerto/electromagnetics/source/plane_wave.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/huerto/electromagnetics/source/beam.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/huerto/electromagnetics/pml/cpml_border.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/huerto/maths/random.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/huerto/maths/functions/core.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/huerto/simulation/task.cpp
)

# Create output directories
set(BIN_DIR ${CMAKE_CURRENT_BINARY_DIR}/bin)
file(MAKE_DIRECTORY ${BIN_DIR})

# Define the dimension-specific flags
set(DIM1_FLAGS -DHUERTO_ONE_DIM)
set(DIM2_FLAGS -DHUERTO_TWO_DIM)
set(DIM3_FLAGS -DHUERTO_THREE_DIM)

# Create the targets for each dimension
foreach(DIMENSION ${DIMENSIONS})
    set(TARGET_NAME mpulse${DIMENSION}d)
    set(DIM_FLAGS DIM${DIMENSION}_FLAGS)
    
    # Create an object library for this dimension's objects
    add_library(${TARGET_NAME}_objects OBJECT ${ALL_SOURCE_FILES})
    
    # Set include directories for this target
    target_include_directories(${TARGET_NAME}_objects PRIVATE ${INCLUDE_DIRS})
    
    # Set compile definitions for this dimension
    target_compile_definitions(${TARGET_NAME}_objects PRIVATE ${${DIM_FLAGS}})
    
    # Create the executable
    add_executable(${TARGET_NAME} $<TARGET_OBJECTS:${TARGET_NAME}_objects>)
    
    # Link libraries
    target_link_libraries(${TARGET_NAME} PRIVATE
        ${HDF5_LIBRARY}
        ${SCHNEK_LIBRARY}
        ${FFTW_LIBRARY}
        ${KOKKOS_CORE_LIBRARY}
        ${KOKKOS_CONTAINERS_LIBRARY}
        ${KOKKOS_SIMD_LIBRARY}
        m
    )
    
    # Set the output directory
    set_target_properties(${TARGET_NAME} PROPERTIES
        RUNTIME_OUTPUT_DIRECTORY ${BIN_DIR}
    )
endforeach()

# Set LD_LIBRARY_PATH for libsz.so issue
# Note: This cannot be done directly in CMake, but we can provide a script
# configure_file(
#     ${CMAKE_CURRENT_SOURCE_DIR}/set_env.sh.in
#     ${CMAKE_CURRENT_BINARY_DIR}/set_env.sh
#     @ONLY
# )

message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")
message(STATUS "Executables will be placed in: ${BIN_DIR}")
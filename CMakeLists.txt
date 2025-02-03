cmake_minimum_required(VERSION 3.10)

project(lattice_processing LANGUAGES Fortran)

set(CMAKE_Fortran_STANDARD 2008)
set(CMAKE_Fortran_STANDARD_REQUIRED ON)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

# Ensure module files are placed correctly
set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/modules)

# Set source directory path
set(SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/source)

# Collect all Fortran source files from the new directory
set(FORTRAN_SOURCES
    # ${SOURCE_DIR}/lattice_processing.f90
    ${SOURCE_DIR}/process_lattice.f90
    # ${SOURCE_DIR}/include/element.f90
    # ${SOURCE_DIR}/include/potential.f90
    ${SOURCE_DIR}/io/atom_io.f90
    ${SOURCE_DIR}/io/input_io.f90
    ${SOURCE_DIR}/io/lattice_io.f90
    ${SOURCE_DIR}/potpars/lmto_lib.f90
    ${SOURCE_DIR}/potpars/potpar_wrapper.f90
    ${SOURCE_DIR}/spheres/defmtr.f90
    ${SOURCE_DIR}/spheres/empty.f90
    ${SOURCE_DIR}/spheres/radii.f90
    ${SOURCE_DIR}/spheres/sphere.f90
    ${SOURCE_DIR}/spheres/utils.f90
    ${SOURCE_DIR}/types/element_info.f90
    ${SOURCE_DIR}/types/system_data.f90
)

# Create the executable
add_executable(process_lattice.x ${FORTRAN_SOURCES})

# Ensure correct module paths
target_include_directories(process_lattice.x PRIVATE ${CMAKE_Fortran_MODULE_DIRECTORY})

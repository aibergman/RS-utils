cmake_minimum_required(VERSION 3.10)

project(process_lattice LANGUAGES Fortran)

# Specify the Fortran standard (e.g., Fortran 2008)
set(CMAKE_Fortran_STANDARD 2008)
set(CMAKE_Fortran_STANDARD_REQUIRED ON)

# Collect all source files from the subdirectories
set(SOURCES
    types/system_data.f90
    io/atom_io.f90
    io/input_io.f90
    io/lattice_io.f90
    spheres/defmtr.f90
    spheres/empty.f90
    spheres/radii.f90
    spheres/sphere.f90
    spheres/utils.f90
    process_lattice.f90
    potpars/potpar_wrapper.f90
    potpars/lmto_lib.f90
)

# Create the executable target
add_executable(process_lattice.x ${SOURCES})

# Specify include directories so that modules can be found.
target_include_directories(process_lattice.x PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/types
    ${CMAKE_CURRENT_SOURCE_DIR}/io
    ${CMAKE_CURRENT_SOURCE_DIR}/spheres
)

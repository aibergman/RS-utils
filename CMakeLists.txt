cmake_minimum_required(VERSION 3.10)
project(rs-util LANGUAGES Fortran)

set(CMAKE_Fortran_STANDARD 2008)
set(CMAKE_Fortran_STANDARD_REQUIRED ON)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/modules)
set(SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/source)

# Enable/disable Python bindings
option(PYTHON_BINDINGS "Enable Python bindings with f90wrap" OFF)

# Collect Fortran source files for the static library
set(LIB_SOURCES
    ${SOURCE_DIR}/rsutil_lib.f90
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

# Create the static library named rsutil_lib
add_library(rsutil_lib STATIC ${LIB_SOURCES})

target_include_directories(rsutil_lib PRIVATE ${CMAKE_Fortran_MODULE_DIRECTORY})

# Collect Fortran source files for the executable
set(EXE_SOURCES
    ${SOURCE_DIR}/process_lattice.f90
)

# Create the Fortran binary named rsutil.x
add_executable(rsutil.x ${EXE_SOURCES})
target_link_libraries(rsutil.x PRIVATE rsutil_lib)
target_include_directories(rsutil.x PRIVATE ${CMAKE_Fortran_MODULE_DIRECTORY})

# Python bindings (if enabled)
if(PYTHON_BINDINGS)
    find_package(Python3 COMPONENTS Interpreter Development REQUIRED)

    # Define the correct output filename from f90wrap
    set(F90WRAP_OUTPUT ${CMAKE_BINARY_DIR}/f90wrap_rsutil_lib.f90)

    # Run f90wrap on rsutil_lib.f90
    add_custom_command(
        OUTPUT ${F90WRAP_OUTPUT}
        COMMAND ${Python3_EXECUTABLE} -m f90wrap ${SOURCE_DIR}/rsutil_lib.f90
        WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
        COMMENT "Generating Fortran-Python bindings with f90wrap..."
        VERBATIM
    )

    # Create a custom target to ensure f90wrap runs before building rsutil
    add_custom_target(f90wrap_generate DEPENDS ${F90WRAP_OUTPUT})

    # Create shared library for Python module rsutil, depending on f90wrap
    add_library(rsutil MODULE ${F90WRAP_OUTPUT})
    add_dependencies(rsutil f90wrap_generate)
    
    target_link_libraries(rsutil PRIVATE rsutil_lib)
    target_include_directories(rsutil PRIVATE ${CMAKE_Fortran_MODULE_DIRECTORY})
endif()

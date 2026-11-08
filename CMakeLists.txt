cmake_minimum_required(VERSION 3.10)
project(rs-util LANGUAGES Fortran)

set(CMAKE_Fortran_STANDARD 2008)
set(CMAKE_Fortran_STANDARD_REQUIRED ON)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/modules)
set(SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/source)

# Enable/disable Python bindings
option(PYTHON_BINDINGS "Enable Python bindings through the ctypes C interface" OFF)

# Collect Fortran source files for the static library
set(LIB_SOURCES
    ${SOURCE_DIR}/types/system_data.f90
    ${SOURCE_DIR}/types/rsutil_config.f90
    ${SOURCE_DIR}/types/element_info.f90
    ${SOURCE_DIR}/rsutil_lib.f90
    ${SOURCE_DIR}/io/atom_io.f90
    ${SOURCE_DIR}/io/input_io.f90
    ${SOURCE_DIR}/io/lattice_io.f90
    ${SOURCE_DIR}/potpars/lmto_lib.f90
    ${SOURCE_DIR}/potpars/potpar_wrapper.f90
    ${SOURCE_DIR}/spheres/defmtr.f90
    ${SOURCE_DIR}/spheres/empty.f90
    ${SOURCE_DIR}/spheres/questaal_sphere_resize.f90
    ${SOURCE_DIR}/spheres/radii.f90
    ${SOURCE_DIR}/spheres/questaal_findes_backend.f90
    ${SOURCE_DIR}/spheres/sphere.f90
    ${SOURCE_DIR}/spheres/utils.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/vendor/questaal_findes/structures.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/vendor/questaal_findes/findes_resize_bridge.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/vendor/questaal_findes/legacy_runtime.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/vendor/questaal_findes/findes.f
)

if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
    set_source_files_properties(${CMAKE_CURRENT_SOURCE_DIR}/vendor/questaal_findes/findes.f PROPERTIES
        COMPILE_OPTIONS "-std=legacy;-ffixed-line-length-none;-fallow-argument-mismatch")
endif()

# Create the static library named rsutil_lib
add_library(rsutil_lib STATIC ${LIB_SOURCES})

target_include_directories(rsutil_lib PRIVATE ${CMAKE_Fortran_MODULE_DIRECTORY})

# Compile the C-compatible entry point once.  Both the standalone executable
# and the ctypes-loaded Python module use the same object, which also avoids
# concurrent Fortran module-file writes during parallel builds.
add_library(rsutil_c_interface OBJECT ${SOURCE_DIR}/rsutil_c_interface.f90)
add_dependencies(rsutil_c_interface rsutil_lib)
target_include_directories(rsutil_c_interface PRIVATE ${CMAKE_Fortran_MODULE_DIRECTORY})

# Collect Fortran source files for the executable
set(EXE_SOURCES
    ${SOURCE_DIR}/process_lattice.f90
)

# Create the Fortran binary named rsutil.x
add_executable(rsutil.x ${EXE_SOURCES} $<TARGET_OBJECTS:rsutil_c_interface>)
target_link_libraries(rsutil.x PRIVATE rsutil_lib)
target_include_directories(rsutil.x PRIVATE ${CMAKE_Fortran_MODULE_DIRECTORY})

# Python bindings (if enabled)
if(PYTHON_BINDINGS)
    find_package(Python3 COMPONENTS Interpreter REQUIRED)

    # Build the package-local shared library directly from the C-compatible
    # entry point used by rsutil/__init__.py through ctypes.
    add_library(rsutil MODULE $<TARGET_OBJECTS:rsutil_c_interface>)
    set_target_properties(rsutil PROPERTIES PREFIX "")
    target_link_libraries(rsutil PRIVATE rsutil_lib)
    target_include_directories(rsutil PRIVATE ${CMAKE_Fortran_MODULE_DIRECTORY})
endif()

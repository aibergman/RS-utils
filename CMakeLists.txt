cmake_minimum_required(VERSION 3.10)

project(lattice_processing LANGUAGES Fortran)

set(CMAKE_Fortran_STANDARD 2008)
set(CMAKE_Fortran_STANDARD_REQUIRED ON)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

# Ensure -fPIC is enabled for shared library
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fPIC")

# Ensure module files are placed correctly
set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/modules)

# Set source directory path
set(SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/source)

# Collect Fortran source files for the library
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

# Create a SHARED library for linking
add_library(rsutil_lib SHARED ${LIB_SOURCES})

# Ensure correct module paths for the library
target_include_directories(rsutil_lib PRIVATE ${CMAKE_Fortran_MODULE_DIRECTORY})

# Install the shared library so Python can find it
set_target_properties(rsutil_lib PROPERTIES PREFIX "" SUFFIX ".so")

# Find Python3 and NumPy
find_package(Python3 REQUIRED COMPONENTS Interpreter Development NumPy)

# Generate f90wrap Python bindings
add_custom_command(
    OUTPUT f90wrap_rsutil_lib.f90
    COMMAND ${Python3_EXECUTABLE} -m numpy.f2py -m rsutil -h f90wrap_rsutil.pyf --overwrite-signature ${SOURCE_DIR}/rsutil_lib.f90
    COMMAND f90wrap -m rsutil_lib ${SOURCE_DIR}/rsutil_lib.f90
    DEPENDS ${SOURCE_DIR}/rsutil_lib.f90
)

# Define f90wrap as a separate build step
add_custom_target(f90wrap_target DEPENDS f90wrap_rsutil_lib.f90)

# Create a Python extension module (Fix: `MODULE`, not `SHARED`)
add_library(rsutil MODULE f90wrap_rsutil_lib.f90 ${LIB_SOURCES})

# Link against rsutil_lib
target_link_libraries(rsutil PRIVATE rsutil_lib)

# Ensure correct module paths for the Python module
target_include_directories(rsutil PRIVATE ${CMAKE_BINARY_DIR}/modules)

# Ensure the Python module is named correctly
set_target_properties(rsutil PROPERTIES PREFIX "" SUFFIX ".so")

# ✅ FIX: Ensure macOS uses the correct flags
if(APPLE)
    set_target_properties(rsutil PROPERTIES LINK_FLAGS "-bundle -undefined dynamic_lookup -fPIC")
else()
    set_target_properties(rsutil PROPERTIES LINK_FLAGS "-shared -fPIC")
endif()

# Ensure f90wrap runs before the final shared library is built
add_dependencies(rsutil f90wrap_target)

# Create the executable and link with rsutil_lib
set(EXE_SOURCES ${SOURCE_DIR}/process_lattice.f90)
add_executable(process_lattice.x ${EXE_SOURCES})
target_link_libraries(process_lattice.x PRIVATE rsutil_lib)

# Ensure correct module paths for the executable
target_include_directories(process_lattice.x PRIVATE ${CMAKE_Fortran_MODULE_DIRECTORY})

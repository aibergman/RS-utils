cmake_minimum_required(VERSION 3.10)

project(lattice_processing LANGUAGES Fortran)

set(CMAKE_Fortran_STANDARD 2008)
set(CMAKE_Fortran_STANDARD_REQUIRED ON)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

# ✅ Ensure -fPIC is enabled for static library linking
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fPIC")

# ✅ Ensure module files are placed correctly
set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/modules)

# ✅ Set source directory path
set(SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/source)

# ✅ Collect Fortran source files for the static library
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

# ✅ Step 1: Create a STATIC library (not SHARED)
add_library(rsutil_lib STATIC ${LIB_SOURCES})

# ✅ Ensure correct module paths for the library
target_include_directories(rsutil_lib PRIVATE ${CMAKE_Fortran_MODULE_DIRECTORY})

# ✅ Install the static library
set_target_properties(rsutil_lib PROPERTIES PREFIX "lib" SUFFIX ".a")

# ✅ Step 2: Create the executable and link with rsutil_lib
set(EXE_SOURCES ${SOURCE_DIR}/process_lattice.f90)
add_executable(process_lattice.x ${EXE_SOURCES})
target_link_libraries(process_lattice.x PRIVATE rsutil_lib)

# ✅ Ensure correct module paths for the executable
target_include_directories(process_lattice.x PRIVATE ${CMAKE_Fortran_MODULE_DIRECTORY})

# 🔹🔹🔹 Python Bindings (Optional) 🔹🔹🔹

# ✅ Add an option to enable Python bindings (default OFF)
#if("$ENV{PIP_BUILD}" STREQUAL "1")
#    set(ENABLE_PYTHON_BINDINGS ON)
#else()
#    option(ENABLE_PYTHON_BINDINGS "Enable Python bindings using f90wrap" OFF)
#endif()
##option(ENABLE_PYTHON_BINDINGS "Enable Python bindings using f90wrap" OFF)
if(DEFINED ENV{FLIT_BUILD})
    message(STATUS "Flit detected: Enabling Python bindings")
    set(ENABLE_PYTHON_BINDINGS ON)
endif()

option(ENABLE_PYTHON_BINDINGS "Enable Python bindings using f90wrap" OFF)

if(ENABLE_PYTHON_BINDINGS)
    message(STATUS "Python bindings enabled")

    # ✅ Python Wrapping: Find Python3 and NumPy
    find_package(Python3 REQUIRED COMPONENTS Interpreter Development NumPy)

    # ✅ Step 3: Ensure f90wrap waits for rsutil_lib to be built
    add_custom_command(
        OUTPUT f90wrap_rsutil_lib.f90
        COMMAND ${Python3_EXECUTABLE} -m numpy.f2py -m rsutil -h f90wrap_rsutil.pyf --overwrite-signature ${SOURCE_DIR}/rsutil_lib.f90
        COMMAND f90wrap -m rsutil_lib ${SOURCE_DIR}/rsutil_lib.f90
        DEPENDS rsutil_lib ${SOURCE_DIR}/rsutil_lib.f90
        WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
    )

    # ✅ Step 4: Define f90wrap as a separate build step
    add_custom_target(f90wrap_target DEPENDS f90wrap_rsutil_lib.f90)

    # ✅ Step 5: Ensure Python module compilation includes the correct module path
    add_custom_command(
        OUTPUT _rsutil.so
        COMMAND ${Python3_EXECUTABLE} -m numpy.f2py --fcompiler=gnu95 -c f90wrap_rsutil_lib.f90 -L. -L${CMAKE_BINARY_DIR} -I${CMAKE_BINARY_DIR}/modules -lrsutil_lib -m _rsutil
        DEPENDS f90wrap_rsutil_lib.f90 rsutil_lib
        WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
    )

    # ✅ Copy the compiled shared library to the `rsutil/` package
    file(GLOB SHARED_LIBRARY "${CMAKE_BINARY_DIR}/_rsutil*.so")

    if(SHARED_LIBRARY)
        file(COPY ${SHARED_LIBRARY} DESTINATION ${CMAKE_SOURCE_DIR}/rsutil)
    endif()

    # ✅ Define Python extension as a build target
    add_custom_target(rsutil_python_extension ALL DEPENDS _rsutil.so)

    # ✅ Ensure f90wrap runs after `rsutil_lib.mod` is created
    add_dependencies(f90wrap_target rsutil_lib)
    add_dependencies(rsutil_python_extension f90wrap_target)
    # ✅ Print the path of the generated shared library to the screen

    # ✅ Ensure compiled `.so` is installed correctly
    install(FILES ${CMAKE_BINARY_DIR}/rsutil/_rsutil${CMAKE_SHARED_LIBRARY_SUFFIX} DESTINATION rsutil)

else()
    message(STATUS "Python bindings disabled (default)")
endif()

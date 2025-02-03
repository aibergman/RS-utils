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

# ✅ Ensure `_rsutil.so` is copied to the correct Python package directory
add_custom_command(
    TARGET rsutil_lib POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy
    ${CMAKE_BINARY_DIR}/_rsutil.so
        ${CMAKE_SOURCE_DIR}/rsutil/_rsutil.so
    COMMENT "Copying _rsutil.so to the Python package directory (rsutil/)"
)

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

        # ✅ Force CMake to Use Virtual Environment Python (Poetry, venv, Conda)
        execute_process(
            COMMAND python3 -c "import sys; print(sys.executable)"
            OUTPUT_VARIABLE VENV_PYTHON
            OUTPUT_STRIP_TRAILING_WHITESPACE
        )
    
        if(NOT VENV_PYTHON)
            message(FATAL_ERROR "No valid Python found in virtual environment.")
        endif()
    
        message(STATUS "Using Python Executable: ${VENV_PYTHON}")

    # ✅ Find Python3 + NumPy
    find_package(Python3 REQUIRED COMPONENTS Interpreter Development.Module NumPy)

    message(STATUS "============================================")
    message(STATUS "Python3_EXECUTABLE         ${Python3_EXECUTABLE}")
    message(STATUS "Python3_SOABI              ${Python3_SOABI}")
    message(STATUS "Python3_INCLUDE_DIRS       ${Python3_INCLUDE_DIRS}")
    message(STATUS "Python3_NumPy_INCLUDE_DIRS ${Python3_NumPy_INCLUDE_DIRS}")
    message(STATUS "============================================")

    if (Python3_FOUND)
        # ✅ Extract NumPy headers
        execute_process(
            COMMAND "${Python3_EXECUTABLE}" -c "import numpy; print(numpy.get_include())"
            OUTPUT_VARIABLE NUMPY_INCLUDE_DIR
            OUTPUT_STRIP_TRAILING_WHITESPACE
        )
        file(TO_CMAKE_PATH ${NUMPY_INCLUDE_DIR} NUMPY_INCLUDE_DIR)

        # ✅ Generate `f90wrap` wrapper files
        add_custom_command(
            OUTPUT f90wrap_rsutil_lib.f90
            COMMAND ${Python3_EXECUTABLE} -m f90wrap -m rsutil_lib ${SOURCE_DIR}/rsutil_lib.f90
            COMMENT "Generating f90wrap Python wrapper"
            DEPENDS ${SOURCE_DIR}/rsutil_lib.f90
        )

        # ✅ Compile Fortran bindings
        add_library(rsutil SHARED f90wrap_rsutil_lib.f90)
        target_link_libraries(rsutil PRIVATE rsutil_lib)
        set_target_properties(rsutil PROPERTIES PREFIX "" OUTPUT_NAME "_rsutil")

                # ✅ Force `.so` instead of `.dylib` (macOS fix)
                set_target_properties(rsutil PROPERTIES
                PREFIX ""                # Remove 'lib' prefix
                SUFFIX ".so"             # Force .so instead of .dylib
                OUTPUT_NAME "_rsutil"    # Ensure correct Python import name
            )

        # ✅ Ensure correct paths
        target_include_directories(rsutil PRIVATE ${CMAKE_Fortran_MODULE_DIRECTORY} ${NUMPY_INCLUDE_DIR})
    else()
        message(FATAL_ERROR "Python environment with NumPy is required for building Python bindings.")
    endif()
else()
    message(STATUS "Python bindings disabled (default)")
endif()



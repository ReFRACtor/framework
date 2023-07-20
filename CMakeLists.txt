# Cmake will complain unless this is before project
cmake_minimum_required(VERSION 3.2)

# This must be called before project() since it sets this variable as empty into the cache
set(CMAKE_BUILD_TYPE Release CACHE STRING "Type of build to create, Release (optimized) by default")

# Define the name of our project
project(refractor)

# Version number
set(REFRACTOR_MAJOR_VERSION "8.0.0")

# Mark that we are using
enable_language(Fortran)

# Force the compiler ID detection early, so it will already know to load GNU.cmake in an effort to detect use of --sysroot under Conda
set(CMAKE_C_COMPILER_ID GNU CACHE STRING "")
set(CMAKE_CXX_COMPILER_ID GNU CACHE STRING "")

# Global setting: Use C++14
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Conda has a separate set of debug flags defined, which is more than just
# adding a -g. If we are doing a debug rather than release build, use those
# flags if we find them.
set(DEBUG_CXXFLAGS "$ENV{DEBUG_CXXFLAGS}" CACHE STRING "Flags to use in Debug")
set(DEBUG_CFLAGS "$ENV{DEBUG_CFLAGS}" CACHE STRING "Flags to use in Debug")
set(DEBUG_FFLAGS "$ENV{DEBUG_FFLAGS}" CACHE STRING "Flags to use in Debug")
if(CMAKE_BUILD_TYPE STREQUAL "Debug")
  if(DEBUG_CXXFLAGS)
    set(CMAKE_CXX_FLAGS "${DEBUG_CXXFLAGS}")
  endif(DEBUG_CXXFLAGS)
  if(DEBUG_CFLAGS)
    set(CMAKE_C_FLAGS "${DEBUG_CFLAGS}")
  endif(DEBUG_CFLAGS)
  if(DEBUG_FFLAGS)
    # Suppress warning about comparing reals. This is often a bad thing
    # to do, but the fortran code we have does this the "right" way, e.g.
    # checking for identical to zero as a flag value.
    set(CMAKE_Fortran_FLAGS "${DEBUG_FFLAGS}  -Wno-compare-reals")
  endif(DEBUG_FFLAGS)
endif(CMAKE_BUILD_TYPE STREQUAL "Debug")

# Add a new build type that adds a few extra optimizations to Release
if(CMAKE_BUILD_TYPE STREQUAL "ReleaseFast")
  set(CMAKE_CXX_FLAGS "$ENV{CXXFLAGS} -Ofast -march=native")
  set(CMAKE_C_FLAGS "$ENV{CFLAGS} -Ofast -march=native")
  # Suppress warning about comparing reals. This is often a bad thing
  # to do, but the fortran code we have does this the "right" way, e.g.
  # checking for identical to zero as a flag value.
  set(CMAKE_Fortran_FLAGS "$ENV{FFLAGS} -Ofast -march=native -Wno-compare-reals -fno-stack-arrays")
endif(CMAKE_BUILD_TYPE STREQUAL "ReleaseFast")

# Anaconda enables the -fopenmp flag for fortran with its standard compiler
# flags. This causes problems with lidort (which was fairly hard to track
# down). -fopenmp implies -frecursive which makes every function in fortran
# recursive. A side effect of this is all variables are declared on the stack
# instead of the heap. lidort uses a large number of variables, which
# completely exhausts the stack. This is particularly hard to track down
# because the effect of the exhausting the stack is a segmentation fault,
# so you don't immediately know the stack is the problem. To avoid this,
# we explicitely turn of openmp even if the flag is initially turn it on.

if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
  add_definitions("-fno-openmp")
endif(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")

# Another annoying thing Anaconda does on the mac is add
# -Wl,-dead_strip_dylibs. Not really clear why this would be imposed on all
# environments, but it is at least currently done. This causes
# "lazy symbol bindings failed" errors, in particular libgsl does lazy
# bindings iwht libgslcblas. It is common to include libraries on the link
# line for reasons other than directly providing symbols, so we really want
# the CMake rule to determine what libraries to include rather than having
# the linker toss libraries out it "thinks" it doesn't need. We strip
# the option out if it gets set by the passed in LDFLAGS from the Anaconda
# environment.
if(CMAKE_EXE_LINKER_FLAGS)
    string(REPLACE "-Wl,-dead_strip_dylibs" "" CMAKE_EXE_LINKER_FLAGS
           ${CMAKE_EXE_LINKER_FLAGS})
endif()

if(CMAKE_MODULE_LINKER_FLAGS)
    string(REPLACE "-Wl,-dead_strip_dylibs" "" CMAKE_MODULE_LINKER_FLAGS
           ${CMAKE_MODULE_LINKER_FLAGS})
endif()

if(CMAKE_SHARED_LINKER_FLAGS)
    string(REPLACE "-Wl,-dead_strip_dylibs" "" CMAKE_SHARED_LINKER_FLAGS
           ${CMAKE_SHARED_LINKER_FLAGS})
endif()

# Disable deprecation warning because luabind uses auto_ptr, enabling C++11 
# causes the warning to be output
add_definitions("-Wno-deprecated-declarations")

# When we are using debug mode, turn on the blitz range checking.
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -DBZ_DEBUG")

# We want to use full paths for all the libraries, otherwise we
# get weird behavior in conda where system libraries get picked out
# instead of the conda library
cmake_policy(SET CMP0060 NEW)

if(POLICY CMP0068)
  cmake_policy(SET CMP0068 NEW)
endif()

# use, i.e. don't skip the full RPATH for the build tree
SET(CMAKE_SKIP_BUILD_RPATH  FALSE)

# By default install into the build directory
if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
    set(CMAKE_INSTALL_PREFIX "${CMAKE_BINARY_DIR}/install" CACHE PATH "..." FORCE)
    message(STATUS "Setting install prefix to: " ${CMAKE_INSTALL_PREFIX})
    message(STATUS "Specify -DCMAKE_INSTALL_PREFIX to cmake to change")
else()
    message(STATUS "Using install prefix: " ${CMAKE_INSTALL_PREFIX})
endif()

SET(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
# For some reason on the Mac we don't automatically figure out that CONDA
# should be on the rpath. So add that. This seems like a bug, the reading of
# CMAKE_INSTALL_RPATH_USE_LINK_PATH seems to indicate it should do this 
# automatically. But it doesn't, so go ahead and set this
if(DEFINED ENV{CONDA_PREFIX})
  list(APPEND CMAKE_INSTALL_RPATH "$ENV{CONDA_PREFIX}/lib")
endif(DEFINED ENV{CONDA_PREFIX})

# add the automatically determined parts of the RPATH
# which point to directories outside the build tree to the install RPATH
SET(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

# Check that required functions exist and define macros to include their use
INCLUDE (CheckFunctionExists)

set(CMAKE_REQUIRED_LIBRARIES ${CMAKE_REQUIRED_LIBRARIES} m)
check_function_exists(feenableexcept HAVE_FEENABLEEXCEPT)
if(HAVE_FEENABLEEXCEPT)
    add_definitions(-DHAVE_FEENABLEEXCEPT)
endif()

# For modules not available in standard cmake modules
list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake/Modules/")

# Expose options that user should supply
set(DEFAULT_ABSCO_DIR "${CMAKE_SOURCE_DIR}/test/unit/data/in/absco")
set(ABSCO_DIR ${DEFAULT_ABSCO_DIR} CACHE PATH "Base path of Absoprtion Coefficient tables are located")
if(ABSCO_DIR STREQUAL DEFAULT_ABSCO_DIR)
    message(WARNING "ABSCO_DIR option not supplied, using internal path with reduced size ABSCO files. Standard unit tests should all pass. Long unit tests using L2_FP_LONG_CHECK will fail without full sized ABSCO files (distributed separately).")
endif()

# Obtain the CM version of the source code
execute_process(COMMAND python ${PROJECT_SOURCE_DIR}/support/refractor/framework/version_util.py ${PROJECT_SOURCE_DIR}
    OUTPUT_VARIABLE REFRACTOR_CM_VERSION OUTPUT_STRIP_TRAILING_WHITESPACE
)

# Place compiled Fortran modules into a common directory
set_property(GLOBAL PROPERTY Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/modules)
include_directories(${CMAKE_BINARY_DIR}/modules)

# Create install directory to silence a warning message about having a
# nonexistent directory on the include line
file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/modules)
 
# Add thirdparty packages to build
add_subdirectory(thirdparty)

# Set up the refractor library
add_subdirectory(lib)

# Add thirdparty dependencies to refractor library
add_dependencies(refractor thirdparty)

# Add the binary tree to the search path for include files so that we will find version.h
# Also include paths set from lib/ dir
include_directories(${PROJECT_BINARY_DIR} ${FP_INCLUDES})
include_directories(SYSTEM ${THIRDPARTY_INCLUDES})

# Include Lua configs and input rules
add_subdirectory(input)

# Add Python binding rules
option(BUILD_PYTHON_BINDING "Build Python SWIG bindings" ON)
option(PIP_DEVELOP_MODE "Install Python in editable mode" OFF)
if(BUILD_PYTHON_BINDING) 
    add_subdirectory(bindings/python)
endif()

# Install support utils
install(DIRECTORY "support/utils/" DESTINATION bin
    USE_SOURCE_PERMISSIONS)

# Install cmake modules
install(DIRECTORY "cmake/Modules" DESTINATION cmake)

# Build documentation for project
add_subdirectory(doc)

# Add unit and full tests
add_subdirectory(test)

# Configure and install enviroment setup script
# Install to base of installation prefix
# Put here at the end so any variables created by subdirectories
# can be included.
configure_file (
    "${PROJECT_SOURCE_DIR}/script/setup_fp_env.sh.in"
    "${PROJECT_BINARY_DIR}/setup_fp_env.sh"
)

install(FILES "${PROJECT_BINARY_DIR}/setup_fp_env.sh" DESTINATION "${CMAKE_INSTALL_PREFIX}")

# Set up cmake project config file install and include locations 
set(INSTALL_CMAKE_DIR "${CMAKE_INSTALL_PREFIX}/cmake")
set(INSTALL_INCLUDE_DIR "${CMAKE_INSTALL_PREFIX}/include")
set(INSTALL_SHARE_DIR "${CMAKE_INSTALL_PREFIX}/share")

# Add refractor (and python if building) targets to the build-tree export set
if(BUILD_PYTHON_BINDING) 
    export(TARGETS refractor _swig_wrap FILE "${PROJECT_BINARY_DIR}/RefractorTargets.cmake")
else()
    export(TARGETS refractor FILE "${PROJECT_BINARY_DIR}/RefractorTargets.cmake")
endif()
export(PACKAGE Refractor)

# Create the RefractorConfig.cmake and RefractorConfigVersion.cmake files
file(RELATIVE_PATH REL_INCLUDE_DIR "${INSTALL_CMAKE_DIR}" "${INSTALL_INCLUDE_DIR}")
file(RELATIVE_PATH REL_SHARE_DIR "${INSTALL_CMAKE_DIR}" "${INSTALL_SHARE_DIR}")
# ... for the build tree
set(CONF_INCLUDE_DIRS "${PROJECT_SOURCE_DIR}" "${PROJECT_BINARY_DIR}")
configure_file(cmake/Config/RefractorConfig.cmake.in "${PROJECT_BINARY_DIR}/RefractorConfig.cmake" @ONLY)
# ... for the install tree
set(CONF_INCLUDE_DIRS "\${REFRACTOR_CMAKE_DIR}/${REL_INCLUDE_DIR}/refractor"
  "\${REFRACTOR_CMAKE_DIR}/${REL_SHARE_DIR}/refractor/swig"
   ${THIRDPARTY_INCLUDES})
configure_file(cmake/Config/RefractorConfig.cmake.in "${PROJECT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/RefractorConfig.cmake" @ONLY)
# ... for both
configure_file(cmake/Config/RefractorConfigVersion.cmake.in  "${PROJECT_BINARY_DIR}/RefractorConfigVersion.cmake" @ONLY)

# Install RefractorConfig.cmake and RefractorConfigVersion.cmake
install(FILES
  "${PROJECT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/RefractorConfig.cmake"
  "${PROJECT_BINARY_DIR}/RefractorConfigVersion.cmake"
  DESTINATION "${INSTALL_CMAKE_DIR}")
# Install the export set for use with the install-tree
install(EXPORT RefractorTargets DESTINATION  "${INSTALL_CMAKE_DIR}")

# Install "binaries" from bin/
install(DIRECTORY ${PROJECT_SOURCE_DIR}/bin/ DESTINATION bin
    USE_SOURCE_PERMISSIONS
    FILES_MATCHING PATTERN "*.py")

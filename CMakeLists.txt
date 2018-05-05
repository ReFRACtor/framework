# This must be called before project() since it sets this variable as empty into the cache
set(CMAKE_BUILD_TYPE Release CACHE STRING "Type of build to create, Release (optimized) by default")

cmake_minimum_required(VERSION 3.2)
enable_language(Fortran)
project(refractor)

# Version number
set(REFRACTOR_MAJOR_VERSION "8.0.0")

# Global setting: Use C++11
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Disable deprecation warning because luabind uses auto_ptr, enabling C++11 
# causes the warning to be output
add_definitions("-Wno-deprecated-declarations")

# Check that required functions exist and define macros to include their use
INCLUDE (CheckFunctionExists)

set(CMAKE_REQUIRED_LIBRARIES ${CMAKE_REQUIRED_LIBRARIES} m)
check_function_exists(feenableexcept HAVE_FEENABLEEXCEPT)
if(HAVE_FEENABLEEXCEPT)
    add_definitions(-DHAVE_FEENABLEEXCEPT)
endif()

# For modules not available in standard cmake modules
list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake/Modules/")

# By default install into the build directory
if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
    set(CMAKE_INSTALL_PREFIX "${CMAKE_BINARY_DIR}/install" CACHE PATH "..." FORCE)
    message(STATUS "Setting install prefix to: " ${CMAKE_INSTALL_PREFIX})
    message(STATUS "Specify -DCMAKE_INSTALL_PREFIX to cmake to change")
else()
    message(STATUS "Using install prefix: " ${CMAKE_INSTALL_PREFIX})
endif()

# Expose options that user should supply
set(ABSCO_DIR $ENV{abscodir} CACHE PATH "Base path of Absoprtion Coefficient tables are located")
if(NOT ABSCO_DIR)
    message(WARNING "ABSCO_DIR option not supplied, tests will fail")
endif()
set(MERRA_DIR $ENV{merradir} CACHE PATH "Base path of MERRA Composite files are located")
if(NOT MERRA_DIR)
    message(WARNING "MERRA_DIR option not supplied, some tests will fail")
endif()

# Obtain the CM version of the source code
execute_process(COMMAND python ${PROJECT_SOURCE_DIR}/support/refractor/version_util.py ${PROJECT_SOURCE_DIR}
    OUTPUT_VARIABLE REFRACTOR_CM_VERSION OUTPUT_STRIP_TRAILING_WHITESPACE
)

# Create version.h file by evaluating cmake variables
configure_file (
    "${PROJECT_SOURCE_DIR}/bin/version.h.in"
    "${PROJECT_BINARY_DIR}/version.h"
)

# Place compiled Fortran modules into a common directory
set_property(GLOBAL PROPERTY Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/modules)
include_directories(${CMAKE_BINARY_DIR}/modules)
 
# Add thirdparty packages to build
add_subdirectory(thirdparty)

# Set up the refractor library
add_subdirectory(lib)

# Add thirdparty dependencies to refractor library
add_dependencies(refractor thirdparty)

# Add the binary tree to the search path for include files so that we will find version.h
# Also include paths set from lib/ dir
include_directories(${PROJECT_BINARY_DIR} ${FP_INCLUDES} ${THIRDPARTY_INCLUDES})

# Create main executable 
add_executable(l2_fp bin/l2_fp.cc)
target_link_libraries(l2_fp refractor ${FP_LIBRARIES} ${THIRDPARTY_LIBRARIES})

# Include Lua configs and input rules
add_subdirectory(input)

# Add Python binding rules
option(BUILD_PYTHON_BINDING "Build Python SWIG bindings" ON)
if(BUILD_PYTHON_BINDING) 
    add_subdirectory(bindings/python)
endif()

# Create installation rule for binary
install(TARGETS l2_fp
        RUNTIME DESTINATION bin
        LIBRARY DESTINATION lib)

# Install support utils
install(DIRECTORY "support/utils/" DESTINATION bin
    USE_SOURCE_PERMISSIONS)

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

# Add refractor (and python if building) targets to the build-tree export set
if(BUILD_PYTHON_BINDING) 
    export(TARGETS refractor _swig_wrap FILE "${PROJECT_BINARY_DIR}/RefractorTargets.cmake")
else()
    export(TARGETS refractor FILE "${PROJECT_BINARY_DIR}/RefractorTargets.cmake")
endif()
export(PACKAGE Refractor)

# Create the RefractorConfig.cmake and RefractorConfigVersion.cmake files
file(RELATIVE_PATH REL_INCLUDE_DIR "${INSTALL_CMAKE_DIR}" "${INSTALL_INCLUDE_DIR}")
MESSAGE("Relative paths and full: ${RELATIVE_PATH} ${REL_INCLUDE_DIR} ${INSTALL_CMAKE_DIR} ${INSTALL_INCLUDE_DIR}")
# ... for the build tree
set(CONF_INCLUDE_DIRS "${PROJECT_SOURCE_DIR}" "${PROJECT_BINARY_DIR}")
configure_file(cmake/Config/RefractorConfig.cmake.in "${PROJECT_BINARY_DIR}/RefractorConfig.cmake" @ONLY)
MESSAGE("CONF_INCLUDE_DIRS for build ${CONF_INCLUDE_DIRS}")
# ... for the install tree
set(CONF_INCLUDE_DIRS "\${REFRACTOR_CMAKE_DIR}/${REL_INCLUDE_DIR}/refractor")
configure_file(cmake/Config/RefractorConfig.cmake.in "${PROJECT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/RefractorConfig.cmake" @ONLY)
MESSAGE("CONF_INCLUDE_DIRS for install ${CONF_INCLUDE_DIRS}")
# ... for both
configure_file(cmake/Config/RefractorConfigVersion.cmake.in  "${PROJECT_BINARY_DIR}/RefractorConfigVersion.cmake" @ONLY)

MESSAGE("\${PROJECT_BINARY_DIR}\${CMAKE_FILES_DIRECTORY} \${INSTALL_CMAKE_DIR} ${PROJECT_BINARY_DIR}${CMAKE_FILES_DIRECTORY} ${INSTALL_CMAKE_DIR}")
# Install RefractorConfig.cmake and RefractorConfigVersion.cmake
install(FILES
  "${PROJECT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/RefractorConfig.cmake"
  "${PROJECT_BINARY_DIR}/RefractorConfigVersion.cmake"
  DESTINATION "${INSTALL_CMAKE_DIR}")
# Install the export set for use with the install-tree
install(EXPORT RefractorTargets DESTINATION  "${INSTALL_CMAKE_DIR}")

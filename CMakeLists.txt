cmake_minimum_required(VERSION 3.2)
enable_language(Fortran)
project(refractor)

# Version number
set (REFRACTOR_MAJOR_VERSION "8.0.0")

# Global setting: Use C++11
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

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

# Obtain the CM version of the source code
execute_process(COMMAND python ${PROJECT_SOURCE_DIR}/support/full_physics/version_util.py ${PROJECT_SOURCE_DIR}
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

# Set up the full_physics library
add_subdirectory(lib)

# Add thirdparty dependencies to full_physics library
add_dependencies(full_physics thirdparty)

# Add the binary tree to the search path for include files so that we will find version.h
# Also include paths set from lib/ dir
include_directories(${PROJECT_BINARY_DIR} ${FP_INCLUDES} ${THIRDPARTY_INCLUDES})

# Create main executable 
add_executable(l2_fp bin/l2_fp.cc)
target_link_libraries(l2_fp full_physics ${FP_LIBRARIES} ${THIRDPARTY_LIBRARIES})

# Create unit test binary, do not include it in the make all target
add_executable(test_all EXCLUDE_FROM_ALL lib/test_all.cc ${FP_TEST_FILES})
target_link_libraries(test_all full_physics ${FP_LIBRARIES} ${THIRDPARTY_LIBRARIES})

# Include Lua configs and input rules
add_subdirectory(input)

# Add Python binding rules
add_subdirectory(bindings/python)

# Create installation rule for binary
install(TARGETS l2_fp
        RUNTIME DESTINATION bin
        LIBRARY DESTINATION lib)

# Install support utils
install(DIRECTORY "support/utils/" DESTINATION bin
    USE_SOURCE_PERMISSIONS)

# Configure and install enviroment setup script
# Install to base of installation prefix
configure_file (
    "${PROJECT_SOURCE_DIR}/script/setup_fp_env.sh.in"
    "${PROJECT_BINARY_DIR}/setup_fp_env.sh"
)

install(FILES "${PROJECT_BINARY_DIR}/setup_fp_env.sh" DESTINATION "${CMAKE_INSTALL_PREFIX}")

# Build documentation for project
add_subdirectory(doc)

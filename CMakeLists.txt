cmake_minimum_required(VERSION 3.2)
enable_language(Fortran)
project(REFRACTOR)

# Global setting: Use C++11
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

######
# Set up version.h file
# Version number
set (REFRACTOR_MAJOR_VERSION "8.0.0")

execute_process(COMMAND python ${PROJECT_SOURCE_DIR}/support/full_physics/version_util.py ${PROJECT_SOURCE_DIR}
    OUTPUT_VARIABLE REFRACTOR_CM_VERSION OUTPUT_STRIP_TRAILING_WHITESPACE
)

# configure a header file to pass some of the CMake settings
# to the source code
configure_file (
    "${PROJECT_SOURCE_DIR}/bin/version.h.in"
    "${PROJECT_BINARY_DIR}/version.h"
)
 
# Create these in parent scope so available for subdirectories to add to them
set(LIBRARY_INCLUDE_DIRS "")
set(LIBRARY_SOURCE_FILES "")

# Call Cmake file in lib/ sub directory
add_subdirectory(lib)

add_library(l2_fp_lib ${LIBRARY_SOURCE_FILES})

# add the binary tree to the search path for include files
# so that we will find version.h
include_directories(${PROJECT_BINARY_DIR} ${LIBRARY_INCLUDE_DIRS})

add_executable(l2_fp bin/l2_fp.cc)
target_link_libraries(l2_fp l2_fp_lib)

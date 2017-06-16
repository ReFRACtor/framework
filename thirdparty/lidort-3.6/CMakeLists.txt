cmake_minimum_required(VERSION 3.2)
enable_language (Fortran)

project(lidort)
set(DIST_VERSION 3.6)

include_directories(${CMAKE_CURRENT_SOURCE_DIR})

file(GLOB_RECURSE LIDORT_SRC "*.F90")

# Configure overridable LIDORT parameters file
if(NOT LIDORT_MAXLAYER)
    set(LIDORT_MAXLAYER 23 CACHE STRING "Maximum number of LIDORT layers")
endif()

if(NOT LIDORT_MAXATMOSWFS)
    set(LIDORT_MAXATMOSWFS 7 CACHE STRING "Maximum number of LIDORT atmospheric weighting functions")
endif()

configure_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/lidort_def/lidort_pars.F90.in"
    "${PROJECT_BINARY_DIR}/lidort_pars.F90"
)
list(APPEND LIDORT_SRC "${PROJECT_BINARY_DIR}/lidort_pars.F90")

add_library(lidort SHARED ${LIDORT_SRC})

set_target_properties(lidort PROPERTIES VERSION ${DIST_VERSION}
                      SOVERSION ${DIST_VERSION})
install(TARGETS lidort
        LIBRARY DESTINATION "${CMAKE_INSTALL_PREFIX}/lib")

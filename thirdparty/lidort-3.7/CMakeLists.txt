cmake_minimum_required(VERSION 3.2)
enable_language (Fortran)

project(lidort)
set(DIST_VERSION 3.7)

include_directories(${CMAKE_CURRENT_SOURCE_DIR})

file(GLOB_RECURSE LIDORT_SRC "*.f90")

# Configure overridable LIDORT parameters file
if(NOT DEFINED LIDORT_MAXLAYER OR "${LIDORT_MAXLAYER}" STREQUAL "")
    set(LIDORT_MAXLAYER 23 CACHE STRING "Maximum number of LIDORT layers" FORCE)
endif()
message(STATUS "Using LIDORT_MAXLAYER value: ${LIDORT_MAXLAYER}")

if(NOT DEFINED LIDORT_MAXATMOSWFS OR "${LIDORT_MAXATMOSWFS}" STREQUAL "")
    set(LIDORT_MAXATMOSWFS 7 CACHE STRING "Maximum number of LIDORT atmospheric weighting functions" FORCE)
endif()
message(STATUS "Using LIDORT_MAXATMOSWFS value: ${LIDORT_MAXATMOSWFS}")

configure_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/lidort_def/lidort_pars.f90.in"
    "${PROJECT_BINARY_DIR}/lidort_pars.f90"
)
list(APPEND LIDORT_SRC "${PROJECT_BINARY_DIR}/lidort_pars.f90")

add_library(lidort SHARED ${LIDORT_SRC})

set_target_properties(lidort PROPERTIES VERSION ${DIST_VERSION}
                      SOVERSION ${DIST_VERSION})
install(TARGETS lidort
        LIBRARY DESTINATION "${CMAKE_INSTALL_PREFIX}/lib")

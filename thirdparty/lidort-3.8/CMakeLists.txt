cmake_minimum_required(VERSION 3.2)
enable_language (Fortran)

project(lidort)
set(DIST_VERSION 3.8)

# Anaconda enables the -fopenmp flag for fortran with its standard compiler
# flags. This causes problems with lidort (which was fairly hard to track
# down). -fopenmp implies -frecursive which makes every function in fortran
# recursive. A side effect of this is all variables are declared on the stack
# instead of the heap. lidort uses a large number of variables, which
# completely exhausts the stack. This is particularly hard to track down
# because the effect of the exhausting the stack is a segmentation fault,
# so you don't immediately know the stack is the problem. To avoid this,
# we explicitely turn of openmp even if the flags initially turn it on.
if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
  add_definitions("-fno-openmp")
endif(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")

include_directories(${CMAKE_CURRENT_SOURCE_DIR})

# Include all Fortran code under this LIDORT directory...
file(GLOB_RECURSE LIDORT_SRC "*.f90")

# Configure overridable LIDORT parameters file
if(NOT DEFINED LIDORT_MAXLAYER OR "${LIDORT_MAXLAYER}" STREQUAL "")
    set(LIDORT_MAXLAYER 70 CACHE STRING "Maximum number of LIDORT layers" FORCE)
endif()
message(STATUS "Using LIDORT_MAXLAYER value: ${LIDORT_MAXLAYER}")

if(NOT DEFINED LIDORT_MAXATMOSWFS OR "${LIDORT_MAXATMOSWFS}" STREQUAL "")
    set(LIDORT_MAXATMOSWFS 7 CACHE STRING "Maximum number of LIDORT atmospheric weighting functions" FORCE)
endif()
message(STATUS "Using LIDORT_MAXATMOSWFS value: ${LIDORT_MAXATMOSWFS}")

# Fill in the values of LIDORT_MAXLAYER and LIDORT_MAXATMOSWFS to be used by LIDORT
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

cmake_minimum_required(VERSION 3.2)
enable_language (Fortran)

project(vlidort)
set(DIST_VERSION 2.8)

# Anaconda enables the -fopenmp flag for fortran with its standard compiler
# flags. This causes problems with vlidort (which was fairly hard to track
# down). -fopenmp implies -frecursive which makes every function in fortran
# recursive. A side effect of this is all variables are declared on the stack
# instead of the heap. vlidort uses a large number of variables, which
# completely exhausts the stack. This is particularly hard to track down
# because the effect of the exhausting the stack is a segmentation fault,
# so you don't immediately know the stack is the problem. To avoid this,
# we explicitely turn of openmp even if the flags initially turn it on.
if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
  add_definitions("-fno-openmp")
endif(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")

include_directories(${CMAKE_CURRENT_SOURCE_DIR})

# Include all Fortran code under this VLIDORT directory...
file(GLOB_RECURSE VLIDORT_SRC "*.f90")

# Configure overridable VLIDORT parameters file
if(NOT DEFINED VLIDORT_MAXLAYER OR "${VLIDORT_MAXLAYER}" STREQUAL "")
    set(VLIDORT_MAXLAYER 70 CACHE STRING "Maximum number of VLIDORT layers" FORCE)
endif()
message(STATUS "Using VLIDORT_MAXLAYER value: ${VLIDORT_MAXLAYER}")

if(NOT DEFINED VLIDORT_MAXATMOSWFS OR "${VLIDORT_MAXATMOSWFS}" STREQUAL "")
    set(VLIDORT_MAXATMOSWFS 7 CACHE STRING "Maximum number of VLIDORT atmospheric weighting functions" FORCE)
endif()
message(STATUS "Using VLIDORT_MAXATMOSWFS value: ${VLIDORT_MAXATMOSWFS}")

# Fill in the values of VLIDORT_MAXLAYER and VLIDORT_MAXATMOSWFS to be used by VLIDORT
configure_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/vlidort_def/vlidort_pars.f90.in"
    "${PROJECT_BINARY_DIR}/vlidort_pars.f90"
)
list(APPEND VLIDORT_SRC "${PROJECT_BINARY_DIR}/vlidort_pars.f90")

add_library(vlidort SHARED ${VLIDORT_SRC})

set_target_properties(vlidort PROPERTIES VERSION ${DIST_VERSION}
                      SOVERSION ${DIST_VERSION})
install(TARGETS vlidort
        LIBRARY DESTINATION "${CMAKE_INSTALL_PREFIX}/lib")

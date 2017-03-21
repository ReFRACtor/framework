cmake_minimum_required(VERSION 3.2)
enable_language (Fortran)

project(lidort)
set(DIST_VERSION 3.6)

include_directories(${CMAKE_CURRENT_SOURCE_DIR})

file(GLOB_RECURSE LIDORT_SRC "*.F90")
MESSAGE("Current directory: ${CMAKE_CURRENT_SOURCE_DIR}")

add_library(lidort SHARED ${LIDORT_SRC})

set_target_properties(lidort PROPERTIES VERSION ${DIST_VERSION}
                      SOVERSION ${DIST_VERSION})
install(TARGETS lidort
        LIBRARY DESTINATION "${CMAKE_INSTALL_PREFIX}/lib")

cmake_minimum_required(VERSION 3.2)
enable_language (Fortran)

project(twostream)

include_directories(${CMAKE_CURRENT_SOURCE_DIR})

file(GLOB TWOSTREAM_SRC "*.F90")
MESSAGE("Current directory: ${CMAKE_CURRENT_SOURCE_DIR}")

add_library(twostream SHARED ${TWOSTREAM_SRC})

install(TARGETS twostream
         LIBRARY DESTINATION "${CMAKE_INSTALL_PREFIX}/lib")

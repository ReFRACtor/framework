cmake_minimum_required(VERSION 3.2)
enable_language (Fortran)

project(first_order)

include_directories(${CMAKE_CURRENT_SOURCE_DIR})

file(GLOB FIRST_ORDER_SRC "*.f90")
MESSAGE("Current directory: ${CMAKE_CURRENT_SOURCE_DIR}")

add_library(first_order SHARED ${FIRST_ORDER_SRC})

install(TARGETS first_order
         LIBRARY DESTINATION "${CMAKE_INSTALL_PREFIX}/lib")

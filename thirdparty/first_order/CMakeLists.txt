cmake_minimum_required(VERSION 3.2)
enable_language (Fortran)

project(first_order)

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

file(GLOB FIRST_ORDER_SRC "*.f90")
MESSAGE("Current directory: ${CMAKE_CURRENT_SOURCE_DIR}")

add_library(first_order SHARED ${FIRST_ORDER_SRC})

install(TARGETS first_order
         LIBRARY DESTINATION "${CMAKE_INSTALL_PREFIX}/lib")

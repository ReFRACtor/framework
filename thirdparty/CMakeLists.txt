cmake_minimum_required(VERSION 3.2)
enable_language(Fortran)

include(ExternalProject)

find_program(CMAKE_LIBTOOL libtool)
find_program(CMAKE_INSTALL install)

# Use these directory for where files get installed
set(THIRDPARTY_INSTALL_PREFIX ${CMAKE_INSTALL_PREFIX})
set(THIRDPARTY_BIN_DIR ${THIRDPARTY_INSTALL_PREFIX}/bin)
set(THIRDPARTY_LIB_DIR ${THIRDPARTY_INSTALL_PREFIX}/lib)
set(THIRDPARTY_INCLUDE_DIR ${THIRDPARTY_INSTALL_PREFIX}/include)

include(blitz.cmake)
include(lua.cmake)
include(luabind.cmake)

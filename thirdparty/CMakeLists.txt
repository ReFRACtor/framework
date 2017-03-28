cmake_minimum_required(VERSION 3.2)
enable_language(Fortran)

include(ExternalProject)

find_program(CMAKE_LIBTOOL libtool)
find_program(CMAKE_INSTALL install)

# Use these directory for where files get installed
if(NOT THIRDPARTY_INSTALL_PREFIX)
    set(THIRDPARTY_INSTALL_PREFIX ${CMAKE_INSTALL_PREFIX})
endif()

set(THIRDPARTY_BIN_DIR ${THIRDPARTY_INSTALL_PREFIX}/bin)
set(THIRDPARTY_LIB_DIR ${THIRDPARTY_INSTALL_PREFIX}/lib)
set(THIRDPARTY_INCLUDE_DIR ${THIRDPARTY_INSTALL_PREFIX}/include)

# Set up where to put Fortran modules
get_property(THIRDPARTY_MODULE_DIR GLOBAL PROPERTY Fortran_MODULE_DIRECTORY)

# Need to check against empty string since the variable is always defined
if("${THIRDPARTY_MODULE_DIR}" STREQUAL "")
    set(THIRDPARTY_MODULE_DIR ${CMAKE_BINARY_DIR}/modules)
endif()

include(blitz/blitz.cmake)
list(APPEND THIRDPARTY_TARGETS blitz)

find_package(Lua)
if(NOT LUA_FOUND)
    include(lua/lua.cmake)
    list(APPEND THIRDPARTY_TARGETS lua)
endif()

find_package(Luabind)
if(NOT LUABIND_FOUND)
    include(luabind/luabind.cmake)
    add_dependencies(luabind lua)
    list(APPEND THIRDPARTY_TARGETS luabind)
endif()

include(2stream/2stream.cmake)
list(APPEND THIRDPARTY_TARGETS twostream)

include(lidort-3.6/lidort.cmake)
list(APPEND THIRDPARTY_TARGETS lidort)

# Custom target that makes it easier to compile all third party independently
add_custom_target(thirdparty DEPENDS ${THIRDPARTY_TARGETS})

# Expose third party variables to parent scope
set(THIRDPARTY_INSTALL_PREFIX ${THIRDPARTY_INSTALL_PREFIX} PARENT_SCOPE)
set(THIRDPARTY_BIN_DIR ${THIRDPARTY_BIN_DIR} PARENT_SCOPE)
set(THIRDPARTY_LIB_DIR ${THIRDPARTY_LIB_DIR} PARENT_SCOPE)
set(THIRDPARTY_INCLUDE_DIR ${THIRDPARTY_INCLUDE_DIR} PARENT_SCOPE)

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
list(APPEND THIRDPARTY_INCLUDES $THIRDPARTY_INCLUDE_DIR)

# Set up where to put Fortran modules
get_property(THIRDPARTY_MODULE_DIR GLOBAL PROPERTY Fortran_MODULE_DIRECTORY)

# Need to check against empty string since the variable is always defined
if("${THIRDPARTY_MODULE_DIR}" STREQUAL "")
    set(THIRDPARTY_MODULE_DIR ${CMAKE_BINARY_DIR}/modules)
endif()

##########################################################
# Find required third party libraries we won't be building

# Set up GSL
# Use hint GSL_ROOT_DIR to specify directory containing GSL
find_package(GSL REQUIRED)
list(APPEND THIRDPARTY_INCLUDES ${GSL_INCLUDES})
list(APPEND THIRDPARTY_LIBRARIES ${GSL_LIBRARIES})

if(GSL_VERSION VERSION_LESS 2.0)
    message(FATAL_ERROR "GSL version > 2.0 required, found ${GSL_VERSION}")
endif()

# Set up HDF5 using the C and C++ bindings
find_package(HDF5 REQUIRED COMPONENTS C CXX)
list(APPEND THIRDPARTY_INCLUDES ${HDF5_C_INCLUDE_DIR})
list(APPEND THIRDPARTY_LIBRARIES ${HDF5_LIBRARIES})

# Set up Boost
find_package(Boost REQUIRED COMPONENTS regex date_time iostreams filesystem system)
list(APPEND THIRDPARTY_INCLUDES ${Boost_INCLUDE_DIRS})
list(APPEND THIRDPARTY_LIBRARIES ${Boost_LIBRARIES})

################################################
# Configure packages we will possibly build here
include(blitz/blitz.cmake)
list(APPEND THIRDPARTY_TARGETS blitz)
list(APPEND THIRDPARTY_LIBRARIES ${THIRDPARTY_LIB_DIR}/libblitz.so)

find_package(Lua)
if(LUA_FOUND)
    list(APPEND THIRDPARTY_INCLUDES ${LUA_INCLUDE_DIR})
    list(APPEND THIRDPARTY_LIBRARIES ${LUA_LIBRARIES})
else()
    include(lua/lua.cmake)
    list(APPEND THIRDPARTY_TARGETS lua)
    list(APPEND THIRDPARTY_LIBRARIES ${THIRDPARTY_LIB_DIR}/liblua.so)
endif()

find_package(Luabind)
if(LUABIND_FOUND)
    list(APPEND THIRDPARTY_INCLUDES ${LUABIND_INCLUDE_DIRS}) 
    list(APPEND THIRDPARTY_LIBRARIES ${LUABIND_LIBRARIES})
else()
    include(luabind/luabind.cmake)
    add_dependencies(luabind lua)
    list(APPEND THIRDPARTY_TARGETS luabind)
endif()

include(2stream/2stream.cmake)
list(APPEND THIRDPARTY_TARGETS twostream)
list(APPEND THIRDPARTY_LIBRARIES ${THIRDPARTY_LIB_DIR}/libtwostream.so)

include(lidort-3.6/lidort.cmake)
list(APPEND THIRDPARTY_TARGETS lidort)
list(APPEND THIRDPARTY_LIBRARIES ${THIRDPARTY_LIB_DIR}/liblidort.so)

# Custom target that makes it easier to compile all third party independently
add_custom_target(thirdparty DEPENDS ${THIRDPARTY_TARGETS})

# Expose third party variables to parent scope
set(THIRDPARTY_INSTALL_PREFIX ${THIRDPARTY_INSTALL_PREFIX} PARENT_SCOPE)
set(THIRDPARTY_BIN_DIR ${THIRDPARTY_BIN_DIR} PARENT_SCOPE)
set(THIRDPARTY_LIB_DIR ${THIRDPARTY_LIB_DIR} PARENT_SCOPE)
set(THIRDPARTY_INCLUDE_DIR ${THIRDPARTY_INCLUDE_DIR} PARENT_SCOPE)
set(THIRDPARTY_INCLUDES ${THIRDPARTY_INCLUDES} PARENT_SCOPE) 
set(THIRDPARTY_LIBRARIES ${THIRDPARTY_LIBRARIES} PARENT_SCOPE)

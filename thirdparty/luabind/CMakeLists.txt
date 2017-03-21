# Copyright (C) 2007-2012 LuaDist.
# Modified by Peter Kapec
# Redistribution and use of this file is allowed according to the terms of the MIT license.
# For details see the COPYRIGHT file distributed with LuaDist.
# Please note that the package source code is licensed under its own license.
# Original build script by Leandro Motta Barros.
#
# Original version from: https://github.com/LuaDist/luabind
# Modified for use by ReFRACtor, not meant be called directly, but copied into 
# luabind extracted archive and run from there

project ( luabind CXX )
cmake_minimum_required ( VERSION 2.8 )
set (DIST_VERSION 0.9.1)

# Find Boost
set ( Boost_USE_STATIC_LIBS OFF )
find_package ( Boost 1.38 REQUIRED )

# Include directories
include_directories ( ${Boost_INCLUDE_DIRS} ${LUA_INCLUDE_DIR} ${CMAKE_CURRENT_SOURCE_DIR} )

# Build the library
file ( GLOB LUABIND_SRC "src/*.cpp" )

add_definitions ( -DLUABIND_BUILDING )

add_library ( luabind SHARED ${LUABIND_SRC} )
target_link_libraries ( luabind ${LUA_LIBRARIES} )

# Install library
set_target_properties ( luabind PROPERTIES VERSION ${DIST_VERSION}
                        SOVERSION ${DIST_VERSION} )
install ( TARGETS luabind
          LIBRARY DESTINATION "${CMAKE_INSTALL_PREFIX}/lib" )

# Install includes from luabind directory
install ( DIRECTORY luabind DESTINATION "${CMAKE_INSTALL_PREFIX}/include" )

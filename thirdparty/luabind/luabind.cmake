set(LUABIND_NAME luabind)

# Using updated and fixed luabind version from:
# https://github.com/Oberon00/luabind
# Commit version used is shown in the filename below
if(NOT DEFINED LUABIND_URL)
  set(LUABIND_URL ${CMAKE_CURRENT_SOURCE_DIR}/luabind/luabind-Oberon00-c0b9359.tar.gz)
endif(NOT DEFINED LUABIND_URL)

# Set up arguments to cmake call
# 1. Tell their cmake that Lua was found and pass the include directory path found by thirdparty cmake
# 2. Turn on shared libary creation
# 3. Don't put the version number in the library filename
# 4. Turn off building of test files, just produce the libary
# 5. Set the install path to the third party install path
set(CMAKE_ARGS -DLUA_FOUND=1 -DLUA_INCLUDE_DIRS=${LUA_INCLUDE_DIR} -DLUABIND_DYNAMIC_LINK=1 -DLUABIND_APPEND_VERSION_SUFFIX=0 -DBUILD_TESTING=0 -DCMAKE_INSTALL_PREFIX=${CMAKE_INSTALL_PREFIX} -DLUA_LIBRARIES=${LUA_LIBRARIES})

# If BOOST_ROOT has been specified for finding boost, then supply this to the third party cmake as well
# Don't want to supply it with an empty location if the variable is not defined
if(BOOST_ROOT)
  list(APPEND CMAKE_ARGS -DBOOST_ROOT=${BOOST_ROOT})
endif(BOOST_ROOT)

# For GCC 7.3, luabind has a segmentation fault when compiled optimized. This
# is probably some bug in luabind that we could fix, but 1) luabind isn't
# really updated anymore and 2) we will move lua out of ReFactor in favor of
# python. So it isn't worth trying to fix. Just turn off optimization for newer
# versions of the compiler. Since this is just used by lua configuration, this
# shouldn't have much of an impact on over all performance
set(LUABIND_CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS})
if(CMAKE_COMPILER_IS_GNUCC AND CMAKE_CXX_COMPILER_VERSION VERSION_GREATER 5.0)
  set(LUABIND_CMAKE_CXX_FLAGS ${LUABIND_CMAKE_CXX_FLAGS} -O0)
endif(CMAKE_COMPILER_IS_GNUCC AND CMAKE_CXX_COMPILER_VERSION VERSION_GREATER 5.0)

# Build and install manually to avoid dependence on bjam
ExternalProject_Add(${LUABIND_NAME}
    URL ${LUABIND_URL}
    CONFIGURE_COMMAND export CXXFLAGS=${CMAKE_CXX_FLAGS} COMMAND export CFLAGS=${CMAKE_C_FLAGS} COMMAND export CXX="${CMAKE_CXX_COMPILER}" COMMAND export CC="${CMAKE_C_COMPILER}" COMMAND cmake ${CMAKE_ARGS} . 
    BUILD_COMMAND make luabind
    BUILD_IN_SOURCE 1
)

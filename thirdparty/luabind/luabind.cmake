set(LUABIND_NAME luabind)

# Using updated and fixed luabind version from:
# https://github.com/Oberon00/luabind
# Commit version used is shown in the filename below
set(LUABIND_URL ${CMAKE_CURRENT_SOURCE_DIR}/luabind/luabind-Oberon00-60e576e.tar.gz)

# Set up arguments to cmake call
# 1. Tell their cmake that Lua was found and pass the include directory path found by thirdparty cmake
# 2. Turn on shared libary creation
# 3. Don't put the version number in the library filename
# 4. Turn off building of test files, just produce the libary
# 5. Set the install path to the third party install path
set(CMAKE_ARGS -DLUA_FOUND=1 -DLUA_INCLUDE_DIRS=${LUA_INCLUDE_DIR} -DLUABIND_DYNAMIC_LINK=1 -DLUABIND_APPEND_VERSION_SUFFIX=0 -DBUILD_TESTING=0 -DCMAKE_INSTALL_PREFIX=${CMAKE_INSTALL_PREFIX})

# If BOOST_ROOT has been specified for finding boost, then supply this to the third party cmake as well
# Don't want to supply it with an empty location if the variable is not defined
if(BOOST_ROOT)
    list(APPEND CMAKE_ARGS -DBOOST_ROOT=${BOOST_ROOT})
endif()

# Build and install manually to avoid dependence on bjam
ExternalProject_Add(${LUABIND_NAME}
    URL ${LUABIND_URL}
    CONFIGURE_COMMAND cmake ${CMAKE_ARGS} . 
    BUILD_COMMAND make luabind
    BUILD_IN_SOURCE 1
)

set(LUABIND_NAME luabind)

set(LUABIND_URL ${CMAKE_CURRENT_SOURCE_DIR}/luabind/luabind-0.9.1.tar.gz)

# This patch was taken from https://gist.github.com/2011636. This fixes a 
# problem with building with boost 1.49 and gcc-4.6.3. This problem is 
# described at http://lists.boost.org/Archives/boost/2012/03/191081.php, but
# basically this just rewrites a #elif as a separate #else + #if. It is unclear
# if this is an actual bug in gcc-4.6.3 or not, but in any case this works
# around this.
set(LUABIND_PATCH1 ${CMAKE_CURRENT_SOURCE_DIR}/luabind/luabind_boost.patch)

# This patch allows luabind to work with Lua 5.2 (it was originally developed
# for Lua 5.1). This originally comes from git://git.colberg.org/luabind.git,
# and is described at http://lua.2524044.n2.nabble.com/Luabind-adapted-to-Lua-5-2-td7582662.html
set(LUABIND_PATCH2 ${CMAKE_CURRENT_SOURCE_DIR}/luabind/luabind_lua5.2.patch)

# This patch lets luabind work with boost >= 1.57. 
# See https://github.com/rpavlik/luabind/pull/23 for details
set(LUABIND_PATCH3 ${CMAKE_CURRENT_SOURCE_DIR}/luabind/luabind_boost_1_57.patch)

# Set up arguments to cmake call
set(CMAKE_ARGS -DLUA_INCLUDE_DIR=${LUA_INCLUDE_DIR} -DCMAKE_INSTALL_PREFIX=${CMAKE_INSTALL_PREFIX})

# Build and install manually to avoid dependence on bjam
ExternalProject_Add(${LUABIND_NAME}
    URL ${LUABIND_URL}
    PATCH_COMMAND patch -p1 < ${LUABIND_PATCH1} &&
        patch -p1 < ${LUABIND_PATCH2} &&
        patch -p1 < ${LUABIND_PATCH3}
        CONFIGURE_COMMAND cp ${CMAKE_CURRENT_SOURCE_DIR}/luabind/CMakeLists.txt . && cmake ${CMAKE_ARGS} . 
    BUILD_IN_SOURCE 1
)

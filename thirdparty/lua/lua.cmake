set(LUA_NAME lua)

set(LUA_URL ${CMAKE_CURRENT_SOURCE_DIR}/lua/lua-5.2.2.tar.gz)

# This patch was taken from gentoo. It add support in the Makefile for building
# a shared library version.
set(LUA_PATCH ${CMAKE_CURRENT_SOURCE_DIR}/lua/lua-make.patch)

# Since these depend on other variables and have a space in them, we
# have to define a new variable so that the values are passed to make below correctly
set(LUA_LIBTOOL "${CMAKE_LIBTOOL} --tag=CC")

# Add a CFLAG necessary for compilation
set(LUA_CFLAGS  "-DLUA_USE_LINUX ${LUA_CFLAGS}")

# By default LUA_LDFLAGS is empty but if its set in our env (e.g. conda) we should use it
set(LUA_LDFLAGS  "$ENV{LDFLAGS}")


ExternalProject_Add(${LUA_NAME}
    URL ${LUA_URL}
    PATCH_COMMAND patch -p1 < ${LUA_PATCH}
    BUILD_IN_SOURCE 1
    CONFIGURE_COMMAND ""
    BUILD_COMMAND cd src && make gentoo_all 
        CXX=${CMAKE_CXX_COMPILER} CC=${CMAKE_C_COMPILER} F77=${CMAKE_Fortran_COMPILER} 
        RPATH=${THIRDPARTY_LIB_DIR}
        LIBTOOL=${LUA_LIBTOOL}
        MYCFLAGS=${LUA_CFLAGS}
        MYLDFLAGS=${LUA_LDFLAGS}
        LUA_LIBS="-lreadline -lncurses"
        LIB_LIBS="-lm -ldl" && cd ..
    INSTALL_COMMAND make gentoo_install INSTALL_TOP=${THIRDPARTY_INSTALL_PREFIX}
        CXX=${CMAKE_CXX_COMPILER} CC=${CMAKE_C_COMPILER} F77=${CMAKE_Fortran_COMPILER} 
        LIBTOOL=${LUA_LIBTOOL}
)

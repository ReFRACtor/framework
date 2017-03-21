set(BLITZ_NAME blitz-0.9)

# Location of tar file to build
set(BLITZ_URL ${CMAKE_CURRENT_SOURCE_DIR}/blitz/${BLITZ_NAME}.tar.gz)

# This patch is described at 
# http://bugs.debian.org/cgi-bin/bugreport.cgi?bug=455661. This adds some
# header files that are needed by gcc >= 4.3. These really should have
# always been there, but the headers in gcc were a bit messy before
set(BLITZ_PATCH ${CMAKE_CURRENT_SOURCE_DIR}/blitz/blitz++.patch)

ExternalProject_Add(${BLITZ_NAME}
    URL ${BLITZ_URL}
    PATCH_COMMAND patch -p1 < ${BLITZ_PATCH}
    CONFIGURE_COMMAND ./configure --prefix ${THIRDPARTY_INSTALL_PREFIX} --enable-shared --enable-shared --disable-doxygen
        CXX=${CMAKE_CXX_COMPILER} CC=${CMAKE_C_COMPILER} F77=${CMAKE_Fortran_COMPILER}
    BUILD_IN_SOURCE 1
    BUILD_COMMAND make 
    INSTALL_COMMAND make install
)

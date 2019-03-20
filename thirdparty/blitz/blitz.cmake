set(BLITZ_NAME blitz)

# Location of tar file to build
set(BLITZ_URL ${CMAKE_CURRENT_SOURCE_DIR}/blitz/blitz-0.10.tar.gz)

# Patch to fix build problem with gcc 7
set(BLITZ_PATCH ${CMAKE_CURRENT_SOURCE_DIR}/blitz/blitz-fix-gcc7.patch)

ExternalProject_Add(${BLITZ_NAME}
    URL ${BLITZ_URL}
    PATCH_COMMAND patch -p1 < ${BLITZ_PATCH}
    CONFIGURE_COMMAND ./configure --prefix ${THIRDPARTY_INSTALL_PREFIX} --enable-shared --enable-shared --disable-doxygen
        CXX=${CMAKE_CXX_COMPILER} CC=${CMAKE_C_COMPILER} F77=${CMAKE_Fortran_COMPILER}
    BUILD_IN_SOURCE 1
    BUILD_COMMAND make 
    INSTALL_COMMAND make install
)

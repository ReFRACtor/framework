ExternalProject_Add(lidort
    URL ${CMAKE_CURRENT_SOURCE_DIR}/lidort-3.8
    CONFIGURE_COMMAND cmake
        -DCMAKE_Fortran_COMPILER=${CMAKE_Fortran_COMPILER}
        # Search our own module directory before any -I in the flags: an
        # installed copy of this library's .mod files (e.g. a conda package
        # under the prefix) would otherwise be found first, since -J is
        # searched after -I, and the build would use stale parameters.
        -DCMAKE_Fortran_FLAGS=-I${THIRDPARTY_MODULE_DIR}\ ${CMAKE_Fortran_FLAGS}
        -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
        -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
        -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
        -DCMAKE_INSTALL_PREFIX=${THIRDPARTY_INSTALL_PREFIX}
        -DCMAKE_Fortran_MODULE_DIRECTORY=${THIRDPARTY_MODULE_DIR}
        -DLIDORT_MAXLAYER=${LIDORT_MAXLAYER}
        -DLIDORT_MAXATMOSWFS=${LIDORT_MAXATMOSWFS}
        -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
        ../lidort
)

set(BLITZ_NAME blitz)

# Location of tar file to build
set(BLITZ_URL ${CMAKE_CURRENT_SOURCE_DIR}/blitz/blitz-39f8859.tar.gz)

# Set up arguments to cmake call
# 1. Use toolchain file to set cmake options that fix compilation issues
# 2. Set the install path to the third party install path
set(CMAKE_ARGS -DCMAKE_TOOLCHAIN_FILE=${CMAKE_CURRENT_SOURCE_DIR}/blitz/toolchain_file.cmake -DCMAKE_INSTALL_PREFIX=${CMAKE_INSTALL_PREFIX} )

# If we are building under conda, we are esentially cross compiling and we want to avoid looking in system paths.
# Further, under conda cmake sets CMAKE_INSTALL_LIBDIR to "lib64" when we expect lib
# We set both: -DCMAKE_FIND_USE_CMAKE_SYSTEM_PATH=FALSE -DCMAKE_INSTALL_LIBDIR=lib to work around.
if(DEFINED ENV{CONDA_PREFIX})
    MESSAGE("Detected source build of blitz under conda. Adding workarounds to CMAKE_ARGS.")
    set(CMAKE_ARGS ${CMAKE_ARGS} -DCMAKE_FIND_USE_CMAKE_SYSTEM_PATH=FALSE -DCMAKE_INSTALL_LIBDIR=lib )
endif()

ExternalProject_Add(${BLITZ_NAME}
    URL ${BLITZ_URL}
    CONFIGURE_COMMAND COMMAND ${CMAKE_COMMAND} ${CMAKE_ARGS} .
    BUILD_IN_SOURCE 1
)

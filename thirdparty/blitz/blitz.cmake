set(BLITZ_NAME blitz)

# Location of tar file to build
set(BLITZ_URL ${CMAKE_CURRENT_SOURCE_DIR}/blitz/blitz-39f8859.tar.gz)

# Set up arguments to cmake call
# 1. Use toolchain file to set cmake options that fix compilation issues
# 2. Set the install path to the third party install path
set(CMAKE_ARGS -DCMAKE_TOOLCHAIN_FILE=${CMAKE_CURRENT_SOURCE_DIR}/blitz/toolchain_file.cmake -DCMAKE_INSTALL_PREFIX=${CMAKE_INSTALL_PREFIX} )

ExternalProject_Add(${BLITZ_NAME}
    URL ${BLITZ_URL}
    CONFIGURE_COMMAND COMMAND ${CMAKE_COMMAND} ${CMAKE_ARGS} .
    BUILD_IN_SOURCE 1
)

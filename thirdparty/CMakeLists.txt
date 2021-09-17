cmake_minimum_required(VERSION 3.2)
enable_language(Fortran)

include(ExternalProject)

find_program(CMAKE_LIBTOOL NAMES glibtool libtool)
find_program(CMAKE_INSTALL NAMES ginstall install)

# Use these directory for where files get installed
if(NOT THIRDPARTY_INSTALL_PREFIX)
    set(THIRDPARTY_INSTALL_PREFIX ${CMAKE_INSTALL_PREFIX})
endif()

set(THIRDPARTY_BIN_DIR ${THIRDPARTY_INSTALL_PREFIX}/bin)
set(THIRDPARTY_LIB_DIR ${THIRDPARTY_INSTALL_PREFIX}/lib)
set(THIRDPARTY_INCLUDE_DIR ${THIRDPARTY_INSTALL_PREFIX}/include)
list(APPEND THIRDPARTY_INCLUDES ${THIRDPARTY_INCLUDE_DIR})

# Create install directory to silence a warning message about having a
# nonexistent directory on the include line
file(MAKE_DIRECTORY ${THIRDPARTY_INCLUDE_DIR})

# Set up where to put Fortran modules
get_property(THIRDPARTY_MODULE_DIR GLOBAL PROPERTY Fortran_MODULE_DIRECTORY)

# Need to check against empty string since the variable is always defined
if("${THIRDPARTY_MODULE_DIR}" STREQUAL "")
    set(THIRDPARTY_MODULE_DIR ${CMAKE_BINARY_DIR}/modules)
endif()

# If using gfortran include the quadmath and gfortran libraries if found
# Not all distributions seem to seperate these out 
if(${CMAKE_Fortran_COMPILER_ID} STREQUAL "GNU")
    find_library(LIB_QUADMATH NAMES quadmath libquadmath)
    if(NOT ${LIB_QUADMATH} MATCHES "-NOTFOUND")
        message(STATUS "Linking with libquadmath: ${LIB_QUADMATH}")
        list(APPEND THIRDPARTY_LIBRARIES ${LIB_QUADMATH})
    endif()

    find_library(LIB_GFORTRAN NAMES gfortran libgfortran)
    if(NOT ${LIB_GFORTRAN} MATCHES "-NOTFOUND")
        message(STATUS "Linking with libgfortran: ${LIB_GFORTRAN}")
        list(APPEND THIRDPARTY_LIBRARIES ${LIB_GFORTRAN})
    endif()
endif()

##########################################################
# Find required third party libraries we won't be building

# Set up GSL
# Use hint GSL_ROOT_DIR to specify directory containing GSL
find_package(GSL REQUIRED)
list(APPEND THIRDPARTY_INCLUDES ${GSL_INCLUDE_DIRS})
list(APPEND THIRDPARTY_LIBRARIES ${GSL_LIBRARIES})

if(GSL_VERSION VERSION_LESS 2.2.1)
    message(FATAL_ERROR "GSL version >= 2.2.1 required, found ${GSL_VERSION}")
endif()

# Set up HDF5 using the C and C++ bindings
find_package(HDF5 REQUIRED COMPONENTS C CXX)
# Different versionso of FindHDF5 may use one or the other of these variable naems, with 
# newer versions saying that HDF5_C_INCLUDE_DIRS deprecates the other
if(HDF5_C_INCLUDE_DIRS)
    list(APPEND THIRDPARTY_INCLUDES ${HDF5_C_INCLUDE_DIRS})
else()
    list(APPEND THIRDPARTY_INCLUDES ${HDF5_C_INCLUDE_DIR})
endif()
list(APPEND THIRDPARTY_LIBRARIES ${HDF5_LIBRARIES})

# Set up Boost

# We want to use the Conda version if it is available, and the user hasn't
# otherwise specified a different version
if(NOT BOOST_ROOT AND DEFINED ENV{CONDA_PREFIX})
  set(BOOST_ROOT $ENV{CONDA_PREFIX})
endif(NOT BOOST_ROOT AND DEFINED ENV{CONDA_PREFIX})
find_package(Boost REQUIRED COMPONENTS regex date_time iostreams filesystem system serialization timer)
list(APPEND THIRDPARTY_INCLUDES ${Boost_INCLUDE_DIRS})
list(APPEND THIRDPARTY_LIBRARIES ${Boost_LIBRARIES})

# boost timer depends on librt, which for some reason don't get added
# by the find package. Look for this, and if we find it add it in.
#find_library(LIB_RT NAMES rt librt HINT ${CONDA_BUILD_SYSROOT}/lib64)
find_library(LIB_RT NAMES rt librt)
if(NOT ${LIB_RT} MATCHES "-NOTFOUND")
  message(STATUS "Linking with librt: ${LIB_RT}")
  list(APPEND THIRDPARTY_LIBRARIES ${LIB_RT})
endif()

# Set up Eigen
find_package(Eigen REQUIRED)
list(APPEND THIRDPARTY_INCLUDES ${EIGEN_INCLUDE_DIRS})

# Set up AER OSS if available
find_package(OSS)
if(OSS_FOUND)
    message(STATUS "AER OSS FM will be built.")
    list(APPEND THIRDPARTY_INCLUDES ${OSS_INCLUDE_DIRS})
    list(APPEND THIRDPARTY_LIBRARIES ${OSS_LIBRARIES})
else()
    message(STATUS "Note: NOT building OSS specific modules")
endif()
# Expose whether we found OSS so we can optionally build our OSS FM
set(OSS_FOUND ${OSS_FOUND} PARENT_SCOPE)

# Set up Google Performance Tools (gperftools) if available
find_package(Gperftools)
if(Gperftools_FOUND)
    message(STATUS "Found Gperftools. Adding to libraries and includes.")
    list(APPEND THIRDPARTY_INCLUDES ${Gperftools_INCLUDE_DIRS})
    list(APPEND THIRDPARTY_LIBRARIES ${Gperftools_LIBRARIES})
else()
    message(STATUS "Did not find Gperftools. Skipping link.")
endif()
# Expose whether we found Gperftools so we can optionally add profiling if present
set(Gperftools_FOUND ${Gperftools_FOUND} PARENT_SCOPE)

################################################
# Configure packages we will possibly build here
find_package(Blitz)
if(BLITZ_FOUND AND (NOT BUILD_BLITZ))
    list(APPEND THIRDPARTY_INCLUDES ${BLITZ_INCLUDE_DIR})
    list(APPEND THIRDPARTY_LIBRARIES ${BLITZ_LIBRARIES})
else()
    message(STATUS "Building blitz from source")
    include(blitz/blitz.cmake)
    list(APPEND THIRDPARTY_TARGETS blitz)
    list(APPEND THIRDPARTY_LIBRARIES ${THIRDPARTY_LIB_DIR}/libblitz${CMAKE_SHARED_LIBRARY_SUFFIX})
    # Override the value put into the CACHE by find_package
    set(BLITZ_INCLUDE_DIR ${THIRDPARTY_INCLUDE_DIR} CACHE PATH "Location of blitz includes" FORCE)
endif()

find_package(Lua 5.2)
if(LUA_FOUND AND (NOT BUILD_LUA))
    list(APPEND THIRDPARTY_INCLUDES ${LUA_INCLUDE_DIR})
    list(APPEND THIRDPARTY_LIBRARIES ${LUA_LIBRARIES})
else()
    message(STATUS "Building Lua from source")
    include(lua/lua.cmake)
    set(LUA_LIBRARIES ${THIRDPARTY_LIB_DIR}/liblua${CMAKE_SHARED_LIBRARY_SUFFIX})
    list(APPEND THIRDPARTY_TARGETS lua)
    list(APPEND THIRDPARTY_LIBRARIES ${THIRDPARTY_LIB_DIR}/liblua${CMAKE_SHARED_LIBRARY_SUFFIX})

    # Override the value put into the CACHE by find_package
    set(LUA_INCLUDE_DIR ${THIRDPARTY_INCLUDE_DIR} CACHE PATH "Location of Lua includes" FORCE)
endif()

find_package(Luabind)
if(LUABIND_FOUND AND (NOT BUILD_LUABIND))
    list(APPEND THIRDPARTY_INCLUDES ${LUABIND_INCLUDE_DIRS}) 
    list(APPEND THIRDPARTY_LIBRARIES ${LUABIND_LIBRARIES})
else()
    message(STATUS "Building Luabind from source")
    include(luabind/luabind.cmake)
    if(BUILD_LUA OR (NOT LUA_FOUND))
        add_dependencies(luabind lua)
    endif()
    list(APPEND THIRDPARTY_TARGETS luabind)
    list(APPEND THIRDPARTY_LIBRARIES ${THIRDPARTY_LIB_DIR}/libluabind${CMAKE_SHARED_LIBRARY_SUFFIX})
endif()

find_package(Twostream)
if(TWOSTREAM_FOUND AND (NOT BUILD_TWOSTREAM))
    list(APPEND THIRDPARTY_LIBRARIES ${TWOSTREAM_LIBRARIES})
else()
    message(STATUS "Building 2stream from source")
    include(2stream/2stream.cmake)
    list(APPEND THIRDPARTY_TARGETS twostream)
    list(APPEND THIRDPARTY_LIBRARIES ${THIRDPARTY_LIB_DIR}/libtwostream${CMAKE_SHARED_LIBRARY_SUFFIX})
endif()


# Duplicates the setting of these properties as done in lidort-*/CMakeCache.txt
# because that cmake file is treated as an external project and not exposed to
# the main cmake enviroment. We want to have these parameter defined at the 
# top level as well, but allow the LIDORT cmake files to still remain independent
#
# Configure overridable LIDORT parameters file
if(NOT DEFINED LIDORT_MAXLAYER OR "${LIDORT_MAXLAYER}" STREQUAL "")
    set(LIDORT_MAXLAYER 70 CACHE STRING "Maximum number of LIDORT layers" FORCE)
endif()

if(NOT DEFINED LIDORT_MAXATMOSWFS OR "${LIDORT_MAXATMOSWFS}" STREQUAL "")
    set(LIDORT_MAXATMOSWFS 7 CACHE STRING "Maximum number of LIDORT atmospheric weighting functions" FORCE)
endif()

# Expose these variables to the main cmake
set(LIDORT_MAXLAYER ${LIDORT_MAXLAYER} PARENT_SCOPE)
set(LIDORT_MAXATMOSWFS ${LIDORT_MAXATMOSWFS} PARENT_SCOPE)

find_package(Lidort)
if(LIDORT_FOUND AND (NOT BUILD_LIDORT))
    list(APPEND THIRDPARTY_LIBRARIES ${LIDORT_LIBRARIES})
else()
    message(STATUS "Building lidort from source")
    include(lidort-3.7/lidort.cmake)
    list(APPEND THIRDPARTY_TARGETS lidort)
    list(APPEND THIRDPARTY_LIBRARIES ${THIRDPARTY_LIB_DIR}/liblidort${CMAKE_SHARED_LIBRARY_SUFFIX})
endif()

find_package(First_order)
if(FIRST_ORDER_FOUND AND (NOT BUILD_FIRST_ORDER))
    list(APPEND THIRDPARTY_LIBRARIES ${FIRST_ORDER_LIBRARIES})
else()
    message(STATUS "Building first_order from source")
    include(first_order/first_order.cmake)
    list(APPEND THIRDPARTY_TARGETS first_order)
    list(APPEND THIRDPARTY_LIBRARIES ${THIRDPARTY_LIB_DIR}/libfirst_order${CMAKE_SHARED_LIBRARY_SUFFIX})
endif()

# Custom target that makes it easier to compile all third party independently
add_custom_target(thirdparty DEPENDS ${THIRDPARTY_TARGETS})

# Expose third party variables to parent scope
set(THIRDPARTY_INSTALL_PREFIX ${THIRDPARTY_INSTALL_PREFIX} PARENT_SCOPE)
set(THIRDPARTY_BIN_DIR ${THIRDPARTY_BIN_DIR} PARENT_SCOPE)
set(THIRDPARTY_LIB_DIR ${THIRDPARTY_LIB_DIR} PARENT_SCOPE)
set(THIRDPARTY_INCLUDE_DIR ${THIRDPARTY_INCLUDE_DIR} PARENT_SCOPE)
set(THIRDPARTY_INCLUDES ${THIRDPARTY_INCLUDES} PARENT_SCOPE) 
set(THIRDPARTY_LIBRARIES ${THIRDPARTY_LIBRARIES} PARENT_SCOPE)

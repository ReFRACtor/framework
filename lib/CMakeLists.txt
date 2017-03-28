# Find required third party libraries

# Set up GSL
# Use hint GSL_ROOT_DIR to specify directory containing GSL
find_package(GSL REQUIRED)
list(APPEND INCLUDES ${GSL_INCLUDES})
list(APPEND LIBS ${GSL_LIBRARIES})

if(GSL_VERSION VERSION_LESS 2.0)
    message(FATAL_ERROR "GSL version > 2.0 required, found ${GSL_VERSION}")
endif()

# Set up HDF5 using the C and C++ bindings
find_package(HDF5 REQUIRED COMPONENTS C CXX)
list(APPEND INCLUDES ${HDF5_C_INCLUDE_DIR})
list(APPEND LIBS ${HDF5_LIBRARIES})

# Set up Boost
find_package(Boost REQUIRED COMPONENTS regex date_time iostreams filesystem system)
list(APPEND INCLUDES ${Boost_INCLUDE_DIRS})
list(APPEND LIBS ${Boost_LIBRARIES})

# Find Lua
find_package(Lua REQUIRED)
list(APPEND INCLUDES ${LUA_INCLUDE_DIR})
list(APPEND LIBS ${LUA_LIBRARIES})
# Tell code that Lua was found, since it is required here, it is always found!
# Needs some more work to make Lua truly optional
add_definitions(-DHAVE_LUA)

# Find LuaBind
find_package(Luabind REQUIRED)
list(APPEND INCLUDES ${LUABIND_INCLUDE_DIRS})
list(APPEND LIBS ${LUABIND_LIBRARIES})

# Find LIDORT
find_library(LIDORT_LIBRARY NAMES lidort HINTS ${THIRDPARTY_LIB_DIR})
list(APPEND LIBS ${LIDORT_LIBRARY})

# Find 2stream
find_library(TWOSTREAM_LIBRARY NAMES twostream HINTS ${THIRDPARTY_LIB_DIR})
list(APPEND LIBS ${TWOSTREAM_LIBRARY})

# Initialize variables used in subdirectories
set(FP_SOURCE_FILES "")

foreach(sub_dir FixedLevel Implementation Interface Lua RegisterOutput Support)
    list(APPEND INCLUDES ${CMAKE_CURRENT_SOURCE_DIR}/${sub_dir})
    add_subdirectory(${sub_dir})
endforeach(sub_dir)

include_directories(${INCLUDES})
add_library(full_physics SHARED ${FP_SOURCE_FILES})
target_link_libraries(full_physics ${LIBS})

# Expose libaries and includes to parent scope
set(FP_INCLUDES ${INCLUDES} PARENT_SCOPE)
set(FP_LIBRARIES ${LIBS} PARENT_SCOPE)

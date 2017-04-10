# Tell code that Lua was found, since it is required here, it is always found!
# Needs some more work to make Lua truly optional
add_definitions(-DHAVE_LUA)

# Initialize variables used in subdirectories
set(FP_SOURCE_FILES "")
set(FP_TEST_FILES "")

foreach(sub_dir FixedLevel Implementation Interface Lua RegisterOutput Support Deprecated)
    list(APPEND FP_INCLUDES ${CMAKE_CURRENT_SOURCE_DIR}/${sub_dir})
    add_subdirectory(${sub_dir})
endforeach(sub_dir)

include_directories(${THIRDPARTY_INCLUDES})
include_directories(${FP_INCLUDES})

add_library(full_physics SHARED ${FP_SOURCE_FILES})
target_link_libraries(full_physics ${THIDRTPARTY_LIBRARIES})

# Expose to the parent scope
set(FP_INCLUDES ${FP_INCLUDES} PARENT_SCOPE)
set(FP_TEST_FILES ${FP_TEST_FILES} PARENT_SCOPE)



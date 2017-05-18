# Tell code that Lua was found, since it is required here, it is always found!
# Needs some more work to make Lua truly optional
add_definitions(-DHAVE_LUA)

# Initialize variables used in subdirectories
set(FP_SOURCE_FILES "")
set(FP_TEST_FILES "")

foreach(sub_dir FixedLevel Implementation Interface Lua RegisterOutput Support Deprecated)
    list(APPEND FP_INCLUDES ${CMAKE_CURRENT_SOURCE_DIR}/${sub_dir})
    add_subdirectory(${sub_dir})

    # Install all include files except excluded patterns such as unit test files
    #
    # The trailing / on the sub_dir makes all files get coppied to the same destination path
    # instead of having the directory (ie Interface/) copied to the destination path
    install(DIRECTORY ${sub_dir}/ DESTINATION include/full_physics
        FILES_MATCHING 
            PATTERN "*.h" 
            PATTERN "*_fixture.h" EXCLUDE
            PATTERN "unit_test*.h" EXCLUDE
            PATTERN "boost" EXCLUDE 
            PATTERN "l_rad_fortran" EXCLUDE
            PATTERN "spectrum_doxygen.h" EXCLUDE)
endforeach(sub_dir)

include_directories(${THIRDPARTY_INCLUDES})
include_directories(${FP_INCLUDES})

add_library(full_physics SHARED ${FP_SOURCE_FILES})
target_link_libraries(full_physics ${THIRDPARTY_LIBRARIES})

# Expose to the parent scope
set(FP_INCLUDES ${FP_INCLUDES} PARENT_SCOPE)
set(FP_TEST_FILES ${FP_TEST_FILES} PARENT_SCOPE)

# Create installation rules
install(TARGETS full_physics
    RUNTIME DESTINATION bin
    LIBRARY DESTINATION lib)

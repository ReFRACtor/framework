# Enable ctest
enable_testing()

# Find some binaries needed by the tests, create a wrapper using these paths
find_program(MAKE_BIN make)
find_program(RUBY_BIN ruby)
find_program(H5LS_BIN h5ls)
get_filename_component(HDF5_BIN ${H5LS_BIN} DIRECTORY)
set(L2_FP_BIN ${PROJECT_BINARY_DIR}/l2_fp)

set(CTEST_WRAP ${CMAKE_CURRENT_BINARY_DIR}/ctest_wrap.sh)
configure_file (
    ctest_wrap.sh.in
    ${CTEST_WRAP}
)

# Create test case fixtures default properties file from CMake variables
configure_file (
    unit/global_fixture_default.cc.in 
    "${CMAKE_CURRENT_BINARY_DIR}/global_fixture_default.cc"
)

set(FP_TEST_FILES ${FP_TEST_FILES} "${CMAKE_CURRENT_BINARY_DIR}/global_fixture_default.cc")

# Create unit test binary, do not include it in the make all target
add_executable(test_all EXCLUDE_FROM_ALL ${FP_TEST_FILES} unit/test_all.cc)
target_link_libraries(test_all full_physics ${FP_LIBRARIES} ${THIRDPARTY_LIBRARIES})

# Define unit tests
foreach(test_file ${FP_TEST_FILES})
    # Gather boost test suite names from test source code
    execute_process(COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/unit/boost_suite_name.sh ${test_file}
        OUTPUT_VARIABLE suite_names)

    # Set up an individual test for each suite
    foreach(suite ${suite_names})
        set(test_name unit/${suite})
        add_test(NAME ${test_name}
            COMMAND ${CMAKE_CURRENT_BINARY_DIR}/test_all --run_test=${suite})
    endforeach()
endforeach()

# Target to run only unit tests
add_custom_target(unit_test COMMAND ${CTEST_WRAP} -R ^unit/ DEPENDS test_all)

# Set up full tests
set(FULL_TEST_RUN_DIR ${CMAKE_CURRENT_BINARY_DIR}/full)
file(MAKE_DIRECTORY ${FULL_TEST_RUN_DIR})


# Get lists of full tests
include(${CMAKE_CURRENT_SOURCE_DIR}/full/test_list.cmake)

foreach(test_name ${FULL_TEST_LIST})
    add_test(NAME full/${test_name}
        COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/full/${test_name}/run_test
        WORKING_DIRECTORY ${FULL_TEST_RUN_DIR})
endforeach()

# Target to run only full tests
add_custom_target(full_test COMMAND ${CTEST_WRAP} -R ^full/ DEPENDS l2_fp)

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

#####
# Create unit test binary, do not include it in the make all target
add_executable(test_all EXCLUDE_FROM_ALL ${FP_TEST_FILES} unit/test_all.cc)
target_link_libraries(test_all full_physics ${FP_LIBRARIES} ${THIRDPARTY_LIBRARIES})

#####
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

#####
# Set up full tests
set(FULL_TEST_RUN_DIR ${CMAKE_CURRENT_BINARY_DIR}/full)
file(MAKE_DIRECTORY ${FULL_TEST_RUN_DIR})

# Get lists of full tests
include(${CMAKE_CURRENT_SOURCE_DIR}/full/test_regular_list.cmake)

foreach(test_name ${FULL_TEST_REGULAR_LIST})
    add_test(NAME full/regular/${test_name}
        COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/full/${test_name}/run_test
        WORKING_DIRECTORY ${FULL_TEST_RUN_DIR})
endforeach()

# Target to run only regular full tests
add_custom_target(full_test_regular COMMAND ${CTEST_WRAP} -R ^full/regular DEPENDS l2_fp)

# Add full test target that depends on regular and Python tests if Python bindings are built
add_custom_target(full_test DEPENDS full_test_regular)

#####
# Python related tests

if(BUILD_PYTHON_BINDING)
    # This target does an install so that Python modules are available to tests
    add_custom_target(python_install
        COMMAND ${CMAKE_COMMAND} --build ${PROJECT_BINARY_DIR} --target install
        DEPENDS python)

    ######
    # Include Python tests if building bindings
    include(${CMAKE_CURRENT_SOURCE_DIR}/full/test_python_list.cmake)

    foreach(test_name ${FULL_TEST_PYTHON_LIST})
        add_test(NAME full/python/${test_name}
            COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/full/${test_name}/run_test
            WORKING_DIRECTORY ${FULL_TEST_RUN_DIR})
    endforeach()

    # Run a Make install so that Python bindings are installed before running these tests
    add_custom_target(full_test_python 
        COMMAND ${CTEST_WRAP} -R ^full/python
        DEPENDS python_install)
    add_dependencies(full_test full_test_python)

    #####
    # Set up nose tests
    set(NOSE_WRAP ${CMAKE_CURRENT_BINARY_DIR}/nose_wrap.sh)
    configure_file(nose_wrap.sh.in ${NOSE_WRAP})

    set(NOSE_TEST_RUN_DIR ${CMAKE_CURRENT_BINARY_DIR}/nose)
    file(MAKE_DIRECTORY ${NOSE_TEST_RUN_DIR})

    add_test(NAME python_nose 
        COMMAND ${NOSE_WRAP} ${PROJECT_SOURCE_DIR}/lib/Python ${PROJECT_SOURCE_DIR}/verification_test
        WORKING_DIRECTORY ${NOSE_TEST_RUN_DIR})

    add_custom_target(nose_test
        COMMAND ${CTEST_WRAP} -R python_nose
        DEPENDS python_install)
endif()

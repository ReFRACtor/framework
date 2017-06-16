# Enable ctest
enable_testing()

# Create unit test binary, do not include it in the make all target
add_executable(test_all EXCLUDE_FROM_ALL unit/test_all.cc ${FP_TEST_FILES})
target_link_libraries(test_all full_physics ${FP_LIBRARIES} ${THIRDPARTY_LIBRARIES})

# Define unit tests
foreach(test_file ${FP_TEST_FILES})
    execute_process(COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/unit/boost_suite_name.sh ${test_file}
        OUTPUT_VARIABLE suite_names)
    foreach(suite ${suite_names})
        set(test_name unit/${suite})
        add_test(NAME ${test_name}
            COMMAND ${CMAKE_CURRENT_BINARY_DIR}/test_all --run_test=${suite})
    endforeach()
endforeach()

# Target to run only unit tests
add_custom_target(unit_test COMMAND ${CMAKE_CTEST_COMMAND} -R ^unit/ DEPENDS test_all)


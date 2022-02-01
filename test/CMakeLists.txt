# Enable ctest
enable_testing()

if(HAVE_LUA)
    # Direct source code to include Lua support
    add_definitions(-DHAVE_LUA)
endif()

# Find some binaries needed by the tests, create a wrapper using these paths
find_program(MAKE_BIN make)
find_program(RUBY_BIN ruby)
find_program(H5LS_BIN h5ls)
get_filename_component(HDF5_BIN ${H5LS_BIN} DIRECTORY)

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

# Define the LIDORT compilation option variables into macros so they can be used by other
# pieces of the code. These are defined and loaded in the thirdparty cmake file
add_definitions(-DLIDORT_MAXLAYER=${LIDORT_MAXLAYER} -DLIDORT_MAXATMOSWFS=${LIDORT_MAXATMOSWFS})

#####
# Create unit test binary, do not include it in the make all target
find_package(Boost REQUIRED COMPONENTS regex date_time iostreams filesystem system serialization timer)
add_executable(test_all EXCLUDE_FROM_ALL ${FP_TEST_FILES} unit/test_all.cc)
target_link_libraries(test_all refractor ${FP_LIBRARIES} ${THIRDPARTY_LIBRARIES})
target_include_directories(test_all PRIVATE ${SWIG_INCLUDE_DIR})
set_target_properties(test_all PROPERTIES INSTALL_RPATH "${PROJECT_BINARY_DIR}/lib;${CMAKE_INSTALL_PREFIX}/lib")
if(DEFINED ENV{CONDA_PREFIX})
  set_property(TARGET test_all APPEND PROPERTY INSTALL_RPATH "$ENV{CONDA_PREFIX}/lib")
endif(DEFINED ENV{CONDA_PREFIX})
set_target_properties(test_all PROPERTIES BUILD_WITH_INSTALL_RPATH TRUE)
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
            COMMAND ${CMAKE_CURRENT_BINARY_DIR}/test_all --log_level=test_suite --run_test=${suite})
    endforeach()
endforeach()

# Target to run only unit tests
add_custom_target(unit_test COMMAND ${CTEST_WRAP} -R ^unit/ DEPENDS test_all)

# Target to run all tests
add_custom_target(all_test COMMAND ${CTEST_WRAP} DEPENDS test_all python_install)

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
add_custom_target(full_test_regular COMMAND ${CTEST_WRAP} -R ^full/regular)

# Add full test target that depends on regular and Python tests if Python bindings are built
add_custom_target(full_test DEPENDS full_test_regular)

#####
# Python related tests

if(BUILD_PYTHON_BINDING)
    # This target does an install so that Python modules are available to tests
    add_custom_target(python_install
        COMMAND ${CMAKE_COMMAND} --build ${PROJECT_BINARY_DIR} --target install
        DEPENDS python)

    set(PYTEST_WRAP ${CMAKE_CURRENT_BINARY_DIR}/pytest_wrap.sh)
    configure_file(pytest_wrap.sh.in ${PYTEST_WRAP})

    set(PYTEST_TEST_RUN_DIR ${CMAKE_CURRENT_BINARY_DIR}/pytest)
    file(MAKE_DIRECTORY ${PYTEST_TEST_RUN_DIR})

    add_test(NAME python_pytest 
        COMMAND ${PYTEST_WRAP} ${PROJECT_SOURCE_DIR}/python/test
        WORKING_DIRECTORY ${PYTEST_TEST_RUN_DIR})

    add_custom_target(python_test
        COMMAND ${CTEST_WRAP} -R python_pytest
        DEPENDS python_install)

endif()

# Tell code that Lua was found, since it is required here, it is always found!
# Needs some more work to make Lua truly optional
add_definitions(-DHAVE_LUA)

# Initialize variables used in subdirectories
set(FP_SOURCE_FILES "")
set(FP_TEST_FILES "")
set(SWIG_INTERFACE_FILES "")
set(SWIG_INCLUDE_DIR "")

# Install all common include files
#
# The top level lib folder has some header and interface files we need to install
file(GLOB LIB_INCLUDE_FILES ${CMAKE_CURRENT_SOURCE_DIR}/*.h)
file(GLOB SWIG_INCLUDE_FILES_ADD ${CMAKE_CURRENT_SOURCE_DIR}/*.i)
list(APPEND SWIG_INCLUDE_FILES ${SWIG_INCLUDE_FILES_ADD})
install(FILES ${LIB_INCLUDE_FILES} DESTINATION include/refractor)

foreach(sub_dir RadiativeTransfer Implementation Interface Lua Support StateMapping Example)
    list(APPEND FP_INCLUDES ${CMAKE_CURRENT_SOURCE_DIR}/${sub_dir})
    add_subdirectory(${sub_dir})

    # Install all include files except excluded patterns such as unit test files
    #
    # The trailing / on the sub_dir makes all files get coppied to the same destination path
    # instead of having the directory (ie Interface/) copied to the destination path
    install(DIRECTORY ${sub_dir}/ DESTINATION include/refractor
        FILES_MATCHING 
            PATTERN "*.h" 
            PATTERN "*.hpp" 
            PATTERN "*_fixture.h" EXCLUDE
            PATTERN "unit_test*.h" EXCLUDE
            PATTERN "l_rad_fortran" EXCLUDE
            PATTERN "spectrum_doxygen.h" EXCLUDE)
    install(DIRECTORY ${sub_dir}/ DESTINATION share/refractor/swig
        FILES_MATCHING 
            PATTERN "*.i")
endforeach(sub_dir)

include_directories(SYSTEM ${THIRDPARTY_INCLUDES})
include_directories(${FP_INCLUDES})

# Set up the swig_rules
add_subdirectory(swig_rules)
include_directories(${SWIG_INCLUDE_DIR})

install(FILES  ${SWIG_EXTRA_HEADER_FILES} DESTINATION include/refractor)
add_library(refractor SHARED ${FP_SOURCE_FILES} ${SWIG_EXTRA_SOURCE_FILES})
target_link_libraries(refractor ${THIRDPARTY_LIBRARIES})
target_compile_definitions(refractor PRIVATE SWIG_MAPPER_NAMESPACE=FullPhysics)
target_compile_definitions(refractor PRIVATE SWIG_HAVE_BOOST_SERIALIZATION)

# Expose to the parent scope
set(FP_INCLUDES ${FP_INCLUDES} PARENT_SCOPE)
set(FP_TEST_FILES ${FP_TEST_FILES} PARENT_SCOPE)
set(SWIG_INTERFACE_FILES ${SWIG_INTERFACE_FILES} PARENT_SCOPE)
set(SWIG_INCLUDE_DIR ${SWIG_INCLUDE_DIR} PARENT_SCOPE)
set(SWIG_INCLUDE_FILES ${SWIG_INCLUDE_FILES} PARENT_SCOPE)

# Install SWIG files
install(FILES ${SWIG_INCLUDE_FILES} DESTINATION share/refractor/swig)

# Create installation rules
install(TARGETS refractor
    EXPORT RefractorTargets
    RUNTIME DESTINATION bin
    LIBRARY DESTINATION lib)

# Find Doxygen and Sphinx so we can build documentation
find_package(Doxygen)

add_custom_target(doc)

if(DOXYGEN_FOUND)
    configure_file(${CMAKE_CURRENT_SOURCE_DIR}/doxygen_api.cfg.in ${CMAKE_CURRENT_BINARY_DIR}/doxygen_api.cfg @ONLY)
    configure_file(${CMAKE_CURRENT_SOURCE_DIR}/doxygen_fortran.cfg.in ${CMAKE_CURRENT_BINARY_DIR}/doxygen_fortran.cfg @ONLY)
    add_custom_command(TARGET doc 
        COMMAND ${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/doxygen_api.cfg
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
        COMMENT "Generating API documentation with Doxygen" VERBATIM)
    add_custom_command(TARGET doc 
        COMMAND ${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/doxygen_fortran.cfg
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
        COMMENT "Generating Fortran documentation with Doxygen" VERBATIM)

    # These are optional so that make install does not fail when the doc target is not run
    install(DIRECTORY "${PROJECT_BINARY_DIR}/doc/api/html/" DESTINATION "doc/api"
        USE_SOURCE_PERMISSIONS OPTIONAL)
    install(DIRECTORY "${PROJECT_BINARY_DIR}/doc/fortran/html/" DESTINATION "doc/fortran"
        USE_SOURCE_PERMISSIONS OPTIONAL)

else()
    message(STATUS "Doxygen not found, will not build API documentation")
endif()

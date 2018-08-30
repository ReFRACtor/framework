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

    # Look for ghp-import to allow importing of Doxygen documentation into gh-pages branch
    find_program(GHP_IMPORT ghp-import)
    if (EXISTS ${GHP_IMPORT})
        add_custom_target(gh-pages)
        add_custom_command(TARGET gh-pages
            COMMAND ${GHP_IMPORT} -n -p "${CMAKE_BINARY_DIR}/doc/api/html"
            COMMENT "Importing API documentation into ghp-pages branch"
            WORKING_DIRECTORY ${CMAKE_SOURCE_DIR})
        add_dependencies(gh-pages doc)
            
    else()
        message(STATUS "ghp-import not found, can not import Doxygen into Github Pages branch")
    endif()

else()
    message(STATUS "Doxygen not found, will not build API documentation")
endif()

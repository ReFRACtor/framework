# Find Python
find_package(PythonLibs 3.0 REQUIRED)
find_package(PythonInterp 3.0 REQUIRED)
include_directories(${PYTHON_INCLUDE_PATH})

# Find SWIG
find_package(SWIG 3.0)

if(${SWIG_FOUND})
    # Include SWIG macros
    include(UseSWIG)

    # Generate swig output into the source tree
    set(CMAKE_SWIG_OUTDIR ${CMAKE_CURRENT_SOURCE_DIR}/swig)

    # Include /lib directory since it also contains some *.i files
    # Inlcude SWIG distribution Python *.i files, including top level dir SWIG_DIR is incorrect
    include_directories(${PROJECT_SOURCE_DIR}/lib ${SWIG_DIR}/python)

    # Find interface files in lib/ directory
    # Specify subdirectories to look into so that the *.i files at the top
    # level are not used since they are meant to be included.
    file(GLOB_RECURSE SWIG_INTERFACE_FILES
        ${PROJECT_SOURCE_DIR}/lib/FixedLevel/*.i
        ${PROJECT_SOURCE_DIR}/lib/Implementation/*.i
        ${PROJECT_SOURCE_DIR}/lib/Interface/*.i
        ${PROJECT_SOURCE_DIR}/lib/RegisterOutput/*.i
        ${PROJECT_SOURCE_DIR}/lib/Support/*.i
        ${PROJECT_SOURCE_DIR}/lib/Deprecated/*.i)

    # Pick and choose Lua directory .i files
    # One of the .i files there is included in the other
    list(APPEND SWIG_INTERFACE_FILES 
        ${PROJECT_SOURCE_DIR}/lib/Lua/lua_state.i
        ${PROJECT_SOURCE_DIR}/lib/Lua/luabind_object.i)

    # Use C++ mode for SWIG
    set_source_files_properties(${SWIG_INTERFACE_FILES} PROPERTIES CPLUSPLUS ON)
    
    # Add the -module option to each *.i file based on it's base name without extention
    foreach(swig_file ${SWIG_INTERFACE_FILES})
        get_filename_component(module_name ${swig_file} NAME_WE)
        set_source_files_properties(${swig_file} PROPERTIES SWIG_FLAGS "-module;${module_name};-cpperraswarn")
    endforeach()

    # Create full_physics SWIG library
    swig_add_module(swig_wrap python ${SWIG_INTERFACE_FILES})
   
    # Link with main library, which is also called full_physics!
    swig_link_libraries(swig_wrap ${PYTHON_LIBRARIES} full_physics)
endif()

# Find where Python distribution libaries should be located by querying distutils
execute_process(COMMAND ${PYTHON_EXECUTABLE} -c "import distutils.sysconfig; print(distutils.sysconfig.get_python_lib(True, False, ''))" OUTPUT_VARIABLE PYTHON_INSTALL_DIR OUTPUT_STRIP_TRAILING_WHITESPACE)

# Install non SWIG Python code into full_physics/
install(DIRECTORY ${PROJECT_SOURCE_DIR}/lib/Python/ DESTINATION ${PYTHON_INSTALL_DIR}/full_physics
    FILES_MATCHING PATTERN "*.py")

# Install SWIG Python code into full_physics_swig/
install(TARGETS _swig_wrap DESTINATION ${PYTHON_INSTALL_DIR}/full_physics_swig)
install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/swig/ DESTINATION ${PYTHON_INSTALL_DIR}/full_physics_swig
    FILES_MATCHING PATTERN "*.py")

install(FILES "${CMAKE_CURRENT_SOURCE_DIR}/__init__.py" DESTINATION ${PYTHON_INSTALL_DIR}/full_physics_swig)

# Expose Python install path to parent scope
set(PYTHON_INSTALL_DIR ${PYTHON_INSTALL_DIR} PARENT_SCOPE)

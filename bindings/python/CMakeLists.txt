# Find Python. Note for conda on the Mac you need to find the interpreter
# before the libs (i.e., the order is important). Otherwise this finds
# the system python for libraries but uses conda for the interpreter, which
# of course causes problems. This goes away if you do this in the other order.
find_package(PythonInterp 3.0 REQUIRED)
find_package(PythonLibs 3.0 REQUIRED)
include_directories(${PYTHON_INCLUDE_PATH})

# Find NumPy
find_package(NumPy REQUIRED)
include_directories(${NUMPY_INCLUDE_DIRS})

# Find SWIG, need at least 3.0.10 because of bugs in
# 3.0.8 related to shared pointers and issues with 3.0.9 importing built code
find_package(SWIG 3.0.10 REQUIRED)

# Check for Doxygen for adding documentation to SWIG files
find_package(Doxygen)

# Generate swig output into the binary directory
set(CMAKE_SWIG_OUTDIR ${CMAKE_CURRENT_BINARY_DIR}/swig)

# Include /lib directory since it also contains some *.i files
# Inlcude SWIG distribution Python *.i files, including top level dir SWIG_DIR is incorrect
include_directories(${PROJECT_SOURCE_DIR}/lib ${SWIG_DIR}/python ${SWIG_INCLUDE_DIR})

# ${SWIG_INTERFACE_FILES} is defined by CMakeLists.txt files in the lib/ directory
# Use SWIG to update generated SWIG files if necessary and to build the library
message(STATUS "Using SWIG directly to build Python bindings")

## Save and reset this after SWIG macros operate so that it does not clean up our SWIG generated files
get_directory_property(orig_addl_make_clean_files ADDITIONAL_MAKE_CLEAN_FILES)

# New policy is to pass SWIG_MODULE_NAME as -module. We don't actually need
# this, but no reason not to turn it on. Without this, you get a warning
# message from CMaker version 3.14 and later
if(POLICY CMP0086)
  cmake_policy(SET CMP0086 NEW)
endif()
if(POLICY CMP0078)
  cmake_policy(SET CMP0078 OLD)
endif()

# Include SWIG macros
include(UseSWIG)

# Process each interface to include documentation if DOXYGEN is available
set(BUILD_USED_SWIG_FILES "")

# Add the -module option to each *.i file based on it's base name without extention
foreach(input_swig_file ${SWIG_INTERFACE_FILES})
    get_filename_component(module_name ${input_swig_file} NAME_WE)
    list(APPEND SWIG_MODULE_NAMES ${module_name})

    if(DOXYGEN_FOUND)
        set(module_out_dir "${CMAKE_CURRENT_BINARY_DIR}/doc/${module_name}")
        set(build_swig_file "${module_out_dir}/${module_name}.i")

        # Create a file that includes both the Doxygen SWIG information and source SWIG interface file
        add_custom_command(OUTPUT "${build_swig_file}"
            COMMAND ${CMAKE_SOURCE_DIR}/swig_rules/interface_file_with_doc.sh ${DOXYGEN_EXECUTABLE} ${module_out_dir} ${input_swig_file} ${build_swig_file}
            COMMENT "Generating documentation for ${module_name}"
            DEPENDS ${input_swig_file})
    else()
        # No Doxygen so use the input_swig_file as the one used in the build
        set(build_swig_file ${input_swig_file})
    endif()

    list(APPEND BUILD_USED_SWIG_FILES ${build_swig_file})
    set_source_files_properties(${build_swig_file} PROPERTIES SWIG_FLAGS "-module;${module_name};-cpperraswarn")
endforeach()

# Use C++ mode for SWIG
set_source_files_properties(${BUILD_USED_SWIG_FILES} PROPERTIES CPLUSPLUS ON)

# Create swig_wrap.cc
file(MAKE_DIRECTORY ${CMAKE_SWIG_OUTDIR})
execute_process(COMMAND ${PYTHON_EXECUTABLE} ${CMAKE_SOURCE_DIR}/swig_rules/create_python_init.py ${CMAKE_SOURCE_DIR}/swig_rules/swig_wrap.tmpl ${CMAKE_SWIG_OUTDIR}/swig_wrap.cc swig refractor ${SWIG_MODULE_NAMES})

# Create refractor SWIG library
if(CMAKE_VERSION VERSION_GREATER 3.7)
    swig_add_library(swig_wrap TYPE SHARED LANGUAGE python SOURCES ${BUILD_USED_SWIG_FILES} ${CMAKE_SWIG_OUTDIR}/swig_wrap.cc)
else()
    # Add module was deprecated at CMake 3.8, only use it for older cmake versions
    swig_add_module(swig_wrap python ${BUILD_USED_SWIG_FILES} ${CMAKE_SWIG_OUTDIR}/swig_wrap.cc)
endif()

# Link with main library, which is also called refractor!

if(APPLE)
  # Compiling with python causes a problem on the Mac using conda.
  # See https://github.com/pytorch/pytorch/commit/73f6715f4725a0723d8171d3131e09ac7abf0666
  # The symbols get resolved when loaded into python, so this doesn't cause
  # a problem.
  swig_link_libraries(swig_wrap refractor)
  set_target_properties(_swig_wrap PROPERTIES LINK_FLAGS "-undefined dynamic_lookup")
else()
  swig_link_libraries(swig_wrap ${PYTHON_LIBRARIES} refractor)
endif(APPLE)
  
# Ensure that only this module has access to the SWIG type information
target_compile_definitions(_swig_wrap PRIVATE SWIG_TYPE_TABLE=refractor)

# Give namespace information for generic swig_rules code
target_compile_definitions(_swig_wrap PRIVATE SWIG_MAPPER_NAMESPACE=FullPhysics)

# Ensure module has .so extension (required by python) on all platforms including OS X (where dylib would otherwise be the case)
set_target_properties(_swig_wrap PROPERTIES SUFFIX ".so")

# Find where Python distribution libaries should be located by querying distutils
execute_process(COMMAND ${PYTHON_EXECUTABLE} -c "import distutils.sysconfig; print(distutils.sysconfig.get_python_lib(True, False, ''))" OUTPUT_VARIABLE PYTHON_INSTALL_DIR OUTPUT_STRIP_TRAILING_WHITESPACE)

# Install non SWIG Python code into refractor/
install(DIRECTORY ${PROJECT_SOURCE_DIR}/lib/Python/ DESTINATION ${PYTHON_INSTALL_DIR}/refractor
  FILES_MATCHING
  PATTERN "*_test.py" EXCLUDE
  PATTERN "test_support/*" EXCLUDE
  PATTERN "test_support" EXCLUDE
  PATTERN "*.py")

# Create version.py file by evaluating cmake variables
configure_file (
    "${CMAKE_CURRENT_SOURCE_DIR}/version.py.in"
    "${CMAKE_CURRENT_BINARY_DIR}/version.py"
)

install(FILES "${CMAKE_CURRENT_BINARY_DIR}/version.py" DESTINATION ${PYTHON_INSTALL_DIR}/refractor)

# Install SWIG Python code into refractor_swig/
install(TARGETS _swig_wrap EXPORT RefractorTargets DESTINATION ${PYTHON_INSTALL_DIR}/refractor_swig)
install(DIRECTORY ${CMAKE_SWIG_OUTDIR}/ DESTINATION ${PYTHON_INSTALL_DIR}/refractor_swig
    FILES_MATCHING PATTERN "*.py")

install(FILES "${CMAKE_SOURCE_DIR}/swig_rules/__init__.py" DESTINATION ${PYTHON_INSTALL_DIR}/refractor_swig)

# On systems with older cmake, it may try to link against lib_swig_wrap.so
#  Add lib_swig_wrap.so as a symlink to _swig_wrap.so to support this
install(CODE "execute_process(COMMAND ${CMAKE_COMMAND} -E create_symlink _swig_wrap.so lib_swig_wrap.so
                              WORKING_DIRECTORY ${CMAKE_INSTALL_PREFIX}/${PYTHON_INSTALL_DIR}/refractor_swig)")

# Expose Python install path to parent scope
set(PYTHON_INSTALL_DIR ${PYTHON_INSTALL_DIR} PARENT_SCOPE)

# Make a custom target to allow doing: make python
add_custom_target(python DEPENDS _swig_wrap)

# Suppress some warning in SWIG code, since we don't generate this ourselves
add_definitions("-Wno-unused-parameter")

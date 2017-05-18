# Files in the current directory to include in the library
set(DIR_SOURCES
    absorber_vmr_fixed_level.cc
    absorber_vmr_fixed_level_output.cc
    absorber_vmr_fixed_level_scaled.cc
    absorber_vmr_fixed_level_scaled_output.cc
    pressure_fixed_level.cc
    pressure_fixed_level_output.cc
    pressure_level_input.cc
    temperature_fixed_level.cc
    temperature_fixed_level_output.cc
)

set(TEST_SOURCES
    absorber_vmr_fixed_level_scaled_test.cc
    absorber_vmr_fixed_level_test.cc
    pressure_fixed_level_output_test.cc
    pressure_fixed_level_test.cc
    pressure_level_input_test.cc
    temperature_fixed_level_output_test.cc
    temperature_fixed_level_test.cc
)

# Full path already included
file(GLOB_RECURSE SWIG_SOURCES *.i)
list(APPEND SWIG_INTERFACE_FILES ${SWIG_SOURCES})

# Add current source directory to source file names
foreach(source_name ${DIR_SOURCES})
    set(FP_SOURCE_FILES ${FP_SOURCE_FILES} ${CMAKE_CURRENT_SOURCE_DIR}/${source_name})
endforeach(source_name)

foreach(source_name ${TEST_SOURCES})
    set(FP_TEST_FILES ${FP_TEST_FILES} ${CMAKE_CURRENT_SOURCE_DIR}/${source_name})
endforeach(source_name)

# Expose the list of source files to the parent cmake file
set(FP_SOURCE_FILES ${FP_SOURCE_FILES} PARENT_SCOPE)
set(FP_TEST_FILES ${FP_TEST_FILES} PARENT_SCOPE)
set(SWIG_INTERFACE_FILES ${SWIG_INTERFACE_FILES} PARENT_SCOPE)

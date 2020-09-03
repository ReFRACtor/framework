# Files in the current directory to include in the library
set(DIR_SOURCES
    state_mapping.cc
    state_mapping_at_indexes.cc
    state_mapping_composite.cc
    state_mapping_gaussian.cc
    state_mapping_interpolate.cc
    state_mapping_linear.cc
    state_mapping_log.cc
    state_mapping_offset.cc
    state_mapping_scale.cc
)

set(TEST_SOURCES
    state_mapping_at_indexes_test.cc
    state_mapping_interpolate_test.cc
    state_mapping_composite_test.cc
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

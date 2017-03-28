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

# Add current source directory to source file names
foreach(source_name ${DIR_SOURCES})
    set(FP_SOURCE_FILES ${FP_SOURCE_FILES} ${CMAKE_CURRENT_SOURCE_DIR}/${source_name})
endforeach(source_name)

# Expose the list of source files to the parent cmake file
set(FP_SOURCE_FILES ${FP_SOURCE_FILES} PARENT_SCOPE)

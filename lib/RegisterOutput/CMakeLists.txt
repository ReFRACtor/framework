# Files in the current directory to include in the library
set(DIR_SOURCES
    absorber_absco_output.cc
    absorber_vmr_met_output.cc
    absorber_vmr_level_output.cc
    absorber_vmr_level_scaled_output.cc
    aerosol_aod_output.cc
    aerosol_consolidated_output.cc
    aerosol_param_output.cc
    altitude_output.cc
    connor_convergence_output.cc
    connor_solver_output.cc
    dispersion_fit_output.cc
    dispersion_polynomial_output.cc
    met_pass_through_output.cc
    empirical_orthogonal_function_output.cc
    error_analysis_output.cc
    fluorescence_effect_output.cc
    forward_model_output.cc
    fts_run_log_output.cc
    gas_vmr_apriori_output.cc
    ground_brdf_output.cc
    ground_coxmunk_output.cc
    ground_coxmunk_plus_lambertian_output.cc
    ground_lambertian_output.cc
    high_res_spectrum_output.cc
    level_1b_output.cc
    max_a_posteriori_output.cc
    oco_forward_model_output.cc
    pressure_output.cc
    radiance_scaling_output.cc
    register_output_base.cc
    source_files_output.cc
    state_vector_output.cc
    stokes_coefficient_fraction_output.cc
    temperature_met_output.cc
    temperature_level_offset_output.cc
    zero_offset_waveform_output.cc
)

set(TEST_SOURCES
    absorber_absco_output_test.cc
    aerosol_aod_output_test.cc
    aerosol_consolidated_output_test.cc
    aerosol_param_output_test.cc
    connor_convergence_output_test.cc
    connor_solver_output_test.cc
    dispersion_polynomial_output_test.cc
    empirical_orthogonal_function_output_test.cc
    error_analysis_output_test.cc
    fluorescence_effect_output_test.cc
    forward_model_output_test.cc
    fts_run_log_output_test.cc
    ground_brdf_output_test.cc
    ground_coxmunk_output_test.cc
    ground_coxmunk_plus_lambertian_output_test.cc
    ground_lambertian_output_test.cc
    high_res_spectrum_output_test.cc
    level_1b_output_test.cc
    oco_forward_model_output_test.cc
    radiance_scaling_output_test.cc
    state_vector_output_test.cc
    stokes_coefficient_fraction_output_test.cc
    temperature_met_output_test.cc
    temperature_level_offset_output_test.cc
    zero_offset_waveform_output_test.cc
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


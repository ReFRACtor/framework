# Files in the current directory to include in the library
set(DIR_SOURCES
    absorber.cc
    absorber_vmr.cc
    absorber_vmr_imp_base.cc
    absorber_vmr_scaled.cc
    aerosol.cc
    aerosol_extinction.cc
    aerosol_extinction_imp_base.cc
    aerosol_property.cc
    aerosol_property_imp_base.cc
    altitude.cc
    composite_initial_guess.cc
    convergence_check.cc
    cost_func.cc
    cost_func_diff.cc
    cost_func_diff_state.cc
    cost_func_state.cc
    cost_function.cc
    cost_minimizer.cc
    dispersion.cc
    forward_model.cc
    forward_model_spectral_grid.cc
    gas_absorption.cc
    generic_object_with_cloud_handling.cc
    ground.cc
    ground_imp_base.cc
    ground_piecewise.cc
    ils.cc
    ils_function.cc
    ils_imp_base.cc
    initial_guess.cc
    instrument.cc
    instrument_correction.cc
    iterative_solver.cc
    iterative_solver_der.cc
    level_1b.cc
    level_1b_info.cc
    level_1b_sample_coefficient.cc
    max_a_posteriori.cc
    max_likelihood.cc
    meteorology.cc
    model_measure.cc
    model_state.cc
    nlls_problem.cc
    nlls_problem_state.cc
    nlls_solver.cc
    noise_model.cc
    observation.cc
    optical_properties.cc
    optical_properties_imp_base.cc
    optical_properties_init_base.cc
    output.cc
    perturbation.cc
    pressure.cc
    pressure_imp_base.cc
    problem_state.cc
    rayleigh.cc
    rayleigh_imp_base.cc
    rt_atmosphere.cc
    sample_grid.cc
    sample_grid_imp_base.cc
    solar_absorption_spectrum.cc
    solar_continuum_spectrum.cc
    solar_doppler_shift.cc
    solar_model.cc
    spectral_window.cc
    spectrum_effect.cc
    spectrum_effect_imp_base.cc
    spectrum_sampling.cc
    stokes_coefficient.cc
    stokes_coefficient_imp_base.cc
    surface_temperature.cc
    temperature.cc
    temperature_imp_base.cc
    temperature_offset.cc
    xsec_table.cc
    xsec_table_imp_base.cc
)

set(TEST_SOURCES
    composite_initial_guess_test.cc
    forward_model_spectral_grid_test.cc
    output_test.cc
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


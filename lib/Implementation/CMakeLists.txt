# Files in the current directory to include in the library
set(DIR_SOURCES
    absco.cc
    absco_aer.cc
    absco_coeff.cc
    absco_hdf.cc
    absorber_absco.cc
    absorber_vmr_level.cc
    absorber_vmr_level_scaled.cc
    absorber_vmr_met.cc
    absorber_xsec.cc
    aerosol_extinction_level.cc
    aerosol_extinction_linear.cc
    aerosol_extinction_log.cc
    aerosol_optical.cc
    aerosol_property_hdf.cc
    aerosol_property_rh_hdf.cc
    aerosol_shape_gaussian.cc
    altitude_hydrostatic.cc
    apply_instrument_units.cc
    atmosphere_legacy.cc
    atmosphere_standard.cc
    bad_sample_noise_model.cc
    bard_ml_problem.cc
    bard_nlls_problem.cc
    brown_nlls_problem.cc
    chisq_convergence.cc
    cloud_3d_effect.cc
    composite_perturbation.cc
    connor_convergence.cc
    connor_cost_function.cc
    connor_iteration_log.cc
    connor_solver.cc
    connor_solver_map.cc
    cost_minimizer_gsl.cc
    dispersion_fit.cc
    dispersion_polynomial.cc
    empirical_orthogonal_function.cc
    error_analysis.cc
    fd_forward_model.cc
    fluorescence_effect.cc
    fm_nlls_problem.cc
    forward_model_with_cloud_handling.cc
    fp_logger.cc
    freudenstein_roth_nlls_problem.cc
    gas_vmr_apriori.cc
    ground_brdf.cc
    ground_coxmunk.cc
    ground_coxmunk_plus_lambertian.cc
    ground_emissivity_piecewise.cc
    ground_emissivity_polynomial.cc
    ground_lambertian.cc
    ground_lambertian_piecewise.cc
    ground_with_cloud_handling.cc
    gsl_lsp.cc
    gsl_mdm.cc
    gsl_sm_lsp.cc
    hdf_constant.cc
    helical_valley_nlls_problem.cc
    ils_gaussian.cc
    ils_grating.cc
    ils_instrument.cc
    ils_fast_apply.cc
    ils_table.cc
    initial_guess_value.cc
    instrument_doppler.cc
    jennrich_sampson_nlls_problem.cc
    level_1b_average.cc
    level_1b_cache.cc
    level_1b_scale_radiance.cc
    log_timing.cc
    max_a_posteriori_sqrt_constraint.cc
    max_a_posteriori_standard.cc
    merra_aerosol.cc
    meyer_ml_problem.cc
    meyer_nlls_problem.cc
    model_measure_bard.cc
    model_measure_meyer.cc
    model_measure_standard.cc
    nlls_max_a_posteriori.cc
    nlls_max_likelihood.cc
    nlls_problem_scaled.cc
    nlls_solver_gsl.cc
    nlls_solver_gsl_lmder.cc
    nlls_solver_gsl_lmsder.cc
    nlls_solver_gsl_sm.cc
    nlls_solver_lm.cc
    nonuniform_spectrum_sampling.cc
    observation_level_1b.cc
    optical_properties_wrt_input.cc
    optical_properties_wrt_rt.cc
    output_hdf.cc
    output_hdf_iteration.cc
    planck.F90
    planck.cc
    powell_nlls_problem.cc
    powell_singular_nlls_problem.cc
    pressure_sigma.cc
    pressure_with_cloud_handling.cc
    radiance_scaling.cc
    radiance_scaling_linear_fit.cc
    radiance_scaling_sv_fit.cc
    radiance_scaling_sv_muses_fit.cc
    rayleigh_bodhaine.cc
    rayleigh_young.cc
    rayleigh_young.f90
    reference_vmr_apriori.cc
    refractive_index.cc
    relative_humidity.cc
    rosenbrock2_nlls_problem.cc
    sample_grid_spectral_domain.cc
    scattering_moment_interpolator.cc
    simple_fixed_spectrum_sampling.cc
    solar_absorption_and_continuum.cc
    solar_absorption_gfit_file.cc
    solar_absorption_gfit_file_f.F90
    solar_absorption_table.cc
    solar_continuum_polynomial.cc
    solar_continuum_table.cc
    solar_doppler_shift_distance_velocity.cc
    solar_doppler_shift_polynomial.cc
    solar_reference_spectrum.cc
    solver_iteration_log.cc
    spectral_window_range.cc
    spectrum_sampling_fixed_spacing.cc
    standard_forward_model.cc
    stokes_coefficient_constant.cc
    stokes_coefficient_fraction.cc
    surface_temperature_direct.cc
    temperature_level.cc
    temperature_level_offset.cc
    temperature_met.cc
    uniform_spectrum_sampling.cc
    xsec_table_simple.cc
    xsec_table_temp_dep.cc
)

set(TEST_SOURCES
    absco_aer_test.cc
    absco_coeff_test.cc
    absco_hdf_test.cc
    absorber_absco_test.cc
    absorber_vmr_level_scaled_test.cc
    absorber_vmr_level_test.cc
    absorber_vmr_met_test.cc
    absorber_xsec_test.cc
    aerosol_extinction_linear_test.cc
    aerosol_extinction_log_test.cc
    aerosol_optical_test.cc
    aerosol_property_hdf_test.cc
    aerosol_property_rh_hdf_test.cc
    aerosol_shape_gaussian_test.cc
    altitude_hydrostatic_test.cc
    atmosphere_fixture.cc
    atmosphere_standard_test.cc
    bad_sample_noise_model_test.cc
    chisq_convergence_test.cc
    cloud_3d_effect_test.cc
    composite_initial_guess_full_guess_test.cc
    composite_perturbation_test.cc
    connor_convergence_test.cc
    connor_cost_function_test.cc
    connor_solver_test.cc
    cost_minimizer_gsl_test.cc
    dispersion_fit_test.cc
    dispersion_polynomial_test.cc
    empirical_orthogonal_function_test.cc
    error_analysis_test.cc
    fluorescence_effect_test.cc
    forward_model_with_cloud_handling_test.cc
    fp_logger_test.cc
    gas_vmr_apriori_test.cc
    ground_brdf_test.cc
    ground_coxmunk_plus_lambertian_test.cc
    ground_coxmunk_test.cc
    ground_emissivity_piecewise_test.cc
    ground_emissivity_polynomial_test.cc
    ground_lambertian_piecewise_test.cc
    ground_lambertian_test.cc
    ground_with_cloud_handling_test.cc
    hdf_constant_test.cc
    ils_gaussian_test.cc
    ils_grating_test.cc
    ils_instrument_test.cc
    ils_fast_apply_test.cc
    ils_table_test.cc
    initial_guess_value_test.cc
    instrument_doppler_test.cc
    level_1b_average_test.cc
    level_1b_cache_test.cc
    level_1b_scale_radiance_test.cc
    merra_aerosol_test.cc
    met_data_fixture.cc
    nlls_max_a_posteriori_test.cc
    nlls_max_likelihood_test.cc
    nlls_solver_gsl_sm_bard_test.cc
    nlls_solver_gsl_sm_brown_test.cc
    nlls_solver_gsl_sm_freudenstein_roth_a_test.cc
    nlls_solver_gsl_sm_freudenstein_roth_b_test.cc
    nlls_solver_gsl_sm_helical_valley_test.cc
    nlls_solver_gsl_sm_jennrich_sampson_test.cc
    nlls_solver_gsl_sm_meyer_test.cc
    nlls_solver_gsl_sm_powell_singular_test.cc
    nlls_solver_gsl_sm_powell_test.cc
    nlls_solver_gsl_sm_rosenbrock2_test.cc
    nlls_solver_gsl_test.cc
    nlls_solver_lm_test.cc
    nonuniform_spectrum_sampling_test.cc
    observation_level_1b_test.cc
    optical_properties_test.cc
    output_hdf_iteration_test.cc
    output_hdf_test.cc
    planck_test.cc
    pressure_sigma_test.cc
    pressure_with_cloud_handling_test.cc
    radiance_scaling_linear_fit_test.cc
    radiance_scaling_sv_fit_test.cc
    radiance_scaling_sv_muses_fit_test.cc
    rayleigh_bodhaine_test.cc
    rayleigh_only_test.cc
    rayleigh_young_test.cc
    reference_vmr_apriori_test.cc
    refractive_index_test.cc
    relative_humidity_test.cc
    sample_grid_spectral_domain_test.cc
    serialized_configuration_fixture.cc
    serialized_configuration_fixture_test.cc
    simple_fixed_spectrum_sampling_test.cc
    solar_absorption_and_continuum_test.cc
    solar_absorption_gfit_file_test.cc
    solar_absorption_table_test.cc
    solar_continuum_polynomial_test.cc
    solar_continuum_table_test.cc
    solar_doppler_shift_distance_velocity_test.cc
    solar_doppler_shift_polynomial_test.cc
    solver_finished_fixture.cc
    spectral_window_range_test.cc
    spectrum_sampling_fixed_spacing_test.cc
    standard_forward_model_test.cc
    stokes_coefficient_constant_test.cc
    stokes_coefficient_fraction_test.cc
    temperature_level_offset_test.cc
    temperature_level_test.cc
    temperature_met_test.cc
    uniform_spectrum_sampling_test.cc
    xsec_table_simple_test.cc
    xsec_table_temp_dep_test.cc
)

# Full path already included
file(GLOB_RECURSE SWIG_SOURCES *.i)
# We remove globbed OSS files below if not building OSS
list(APPEND SWIG_INTERFACE_FILES ${SWIG_SOURCES})

if(OSS_FOUND)
    message(STATUS "Adding AER OSS FM to sources")
    list(APPEND DIR_SOURCES oss_forward_model.cc)
    list(APPEND DIR_SOURCES oss_interface.cc)
    list(APPEND DIR_SOURCES oss_retrieval_flags.cc)
    list(APPEND TEST_SOURCES oss_configuration_fixture.cc)
    list(APPEND TEST_SOURCES oss_forward_model_test.cc)
    list(APPEND TEST_SOURCES oss_interface_test.cc)
else()
    list(REMOVE_ITEM SWIG_INTERFACE_FILES ${CMAKE_CURRENT_SOURCE_DIR}/oss_forward_model.i)
    list(REMOVE_ITEM SWIG_INTERFACE_FILES ${CMAKE_CURRENT_SOURCE_DIR}/oss_interface.i)
    list(REMOVE_ITEM SWIG_INTERFACE_FILES ${CMAKE_CURRENT_SOURCE_DIR}/oss_retrieval_flags.i)
    list(REMOVE_ITEM SWIG_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/oss_forward_model.i)
    list(REMOVE_ITEM SWIG_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/oss_interface.i)
    list(REMOVE_ITEM SWIG_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/oss_retrieval_flags.i)
endif()

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


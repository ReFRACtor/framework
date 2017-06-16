# Files in the current directory to include in the library
set(DIR_SOURCES
    acos_ecmwf.cc
    acos_sounding_id.cc
    array_with_unit.cc
    boost/backtrace.cpp
    constant.cc
    constants.F90
    default_constant.cc
    double_with_unit.cc
    environment_substitute.cc
    fp_gsl_integrate.cc
    fp_gsl_matrix.cc
    fp_time.cc
    fstream_compress.cc
    fts_run_log.cc
    hdf_file.cc
    hdf_file_generating.cc
    hdf_sounding_id.cc
    heritage_file.cc
    heritage_matrix_write.cc
    ifstream_cs.cc
    linear_algebra.cc
    logger.cc
    logger_f.F90
    oco_met_file.cc
    oco_sim_met_ecmwf.cc
    oco_sounding_id.cc
    old_constant.cc
    polynomial_eval.cc
    rayleigh_greek_moment.cc
    rf_gauleg.cc
    spectral_bound.cc
    spectral_bound.h
    spectral_domain.cc
    spectral_range.cc
    spectrum.cc
    state_vector.cc
    sub_state_vector_proxy.cc
    unit.cc
    uq_ecmwf.cc
    uq_sounding_id.cc
)

set(TEST_SOURCES
    acos_sounding_id_test.cc
    array_ad_cache_test.cc
    array_ad_test.cc
    auto_derivative_test.cc
    bin_map_test.cc
    closest_point_test.cc
    ecmwf_test.cc
    environment_substitute_test.cc
    fp_exception_test.cc
    fp_gsl_integrate_test.cc
    fp_gsl_matrix_test.cc
    fp_time_test.cc
    fstream_compress_test.cc
    fts_run_log_test.cc
    global_fixture.cc
    hdf_file_test.cc
    heritage_file_test.cc
    heritage_matrix_write_test.cc
    ifstream_cs_test.cc
    linear_algebra_test.cc
    linear_interpolate_test.cc
    log_interpolate_test.cc
    logger_test.cc
    observer_test.cc
    oco_sim_met_ecmwf_test.cc
    oco_sounding_id_test.cc
    ostream_pad_test.cc
    polynomial_eval_test.cc
    printable_test.cc
    rayleigh_greek_moment_test.cc
    spectral_bound_test.cc
    spectral_domain_test.cc
    state_vector_test.cc
    unit_test.cc
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

# Files in the current directory to include in the library
set(DIR_SOURCES
    array_ad.cc
    array_ad_with_unit.cc
    array_with_unit.cc
    auto_derivative.cc
    auto_derivative_with_unit.cc
    boost/backtrace.cpp
    blitz_array_serialize_support.cc
    boost_rational_serialize_support.cc
    generic_object.cc
    generic_object_map.cc
    constant.cc
    constants.F90
    default_constant.cc
    double_with_unit.cc
    mapping.cc
    mapping_gaussian.cc
    mapping_linear.cc
    mapping_log.cc
    mapping_offset.cc
    mapping_scale.cc
    environment_substitute.cc
    fp_gsl_integrate.cc
    fp_gsl_matrix.cc
    fp_time.cc
    fstream_compress.cc
    hdf_file.cc
    hdf_file_generating.cc
    ifstream_cs.cc
    linear_algebra.cc
    logger.cc
    logger_f.F90
    mapping_linear.h
    mapping_log.h
    mapping_scale.h
    mapping_gaussian.h
    old_constant.cc
    polynomial_eval.cc
    rayleigh_greek_moment.cc
    rf_gauleg.cc
    spectral_bound.cc
    spectral_bound.h
    spectral_domain.cc
    spectral_range.cc
    named_spectrum.cc
    spectrum.cc
    stacked_radiance_mixin.cc
    state_vector.cc
    state_vector_observer.cc
    sub_state_vector_observer.cc
    sub_state_vector_proxy.cc
    unit.cc
)

set(TEST_SOURCES
    angle_util_test.cc
    array_ad_cache_test.cc
    array_ad_test.cc
    auto_derivative_test.cc
    bin_map_test.cc
    closest_point_test.cc
    environment_substitute_test.cc
    fp_exception_test.cc
    fp_gsl_integrate_test.cc
    fp_gsl_matrix_test.cc
    fp_time_test.cc
    fstream_compress_test.cc
    global_fixture.cc
    hdf_file_test.cc
    ifstream_cs_test.cc
    linear_algebra_test.cc
    linear_interpolate_test.cc
    log_interpolate_test.cc
    logger_test.cc
    observer_test.cc
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

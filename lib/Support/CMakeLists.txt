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
    ecmwf.cc
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
    oco_ecmwf.cc
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

# Add current source directory to source file names
foreach(source_name ${DIR_SOURCES})
    set(FP_SOURCE_FILES ${FP_SOURCE_FILES} ${CMAKE_CURRENT_SOURCE_DIR}/${source_name})
endforeach(source_name)

# Expose the list of source files to the parent cmake file
set(FP_SOURCE_FILES ${FP_SOURCE_FILES} PARENT_SCOPE)

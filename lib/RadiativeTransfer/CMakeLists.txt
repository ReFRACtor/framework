# Files in the current directory to include in the library
set(DIR_SOURCES
    l_rad_fortran/calc_geom_first_enh.F90
    l_rad_fortran/calc_geom_first_reg.F90
    l_rad_fortran/calc_geom_second.F90
    l_rad_fortran/l_calc_phase_first.F90
    l_rad_fortran/l_rad_first_enh.F90
    l_rad_fortran/l_rad_first_reg.F90
    l_rad_fortran/l_rad_second.F90
    l_rad_fortran/l_surface.F90
    l_rad_fortran/l_surface_fourier.F90
    l_rad_fortran/geometry_pool.F90
    l_rad_fortran/l_rad_first_pp.F90
    l_rad_fortran/l_sscorr_nadir.F90
    l_rad_fortran/l_sscorr_outgoing.F90
    l_rad_fortran/l_sscorr_pp.F90
    pca_fortran/pca_auxiliaries.f90
    pca_fortran/pca_eigensolvers.f90
    pca_fortran/UVVSWIR_Eofpc_CreateBin.f90
    chapman_boa.cc
    chapman_boa_rt.cc
    first_order_interface.F90
    first_order_driver.cc
    first_order_rt.cc
    hres_wrapper.cc
    l_rad_driver.cc
    l_rad_driver_f.F90
    l_rad_rt.cc
    lidort_driver.cc
    lidort_interface_masters.F90
    lidort_interface_masters_io.F90
    lidort_interface_types.F90
    lidort_interface_types_io.F90
    lidort_rt.cc
    lsi_rt.cc
    pca_binning.cc
    pca_eigensolver.cc
    pca_optical_properties.cc
    radiant2.8_direct.F90
    radiant_driver.cc
    radiative_transfer.cc
    radiative_transfer_fixed_stokes_coefficient.cc
    radiative_transfer_single_wn.cc
    spurr_driver.cc
    spurr_rt.cc
    twostream_driver.cc
    twostream_interface.F90
    twostream_rt.cc
)

set(TEST_SOURCES
    chapman_boa_rt_test.cc
    chapman_boa_test.cc
    first_order_driver_test.cc
    first_order_rt_test.cc
    l_rad_driver_test.cc
    l_rad_rt_test.cc
    lidort_driver_test.cc
    lidort_fixture.cc
    lidort_interface_masters_test.cc
    lidort_interface_types_test.cc
    lidort_rt_test.cc
    lsi_rt_test.cc
    pca_binning_test.cc
    pca_eigensolver_test.cc
    twostream_driver_test.cc
    twostream_interface_test.cc
    twostream_rt_test.cc
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


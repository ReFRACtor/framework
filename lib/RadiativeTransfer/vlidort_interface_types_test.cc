#include <blitz/array.h>
#include <cmath>
#include "unit_test_support.h"
#include "vlidort_interface_types.h"


using namespace FullPhysics;
using namespace blitz;

/* This file was auto-generated */

BOOST_FIXTURE_TEST_SUITE(vlidort_interface_types, GlobalFixture)


BOOST_AUTO_TEST_CASE(vlidort_pars)
{
  // Test obtaining instance
  VLidort_Pars tst_obj = VLidort_Pars::instance();
  // Test that the value obtained from fortran is as expected
  BOOST_CHECK_EQUAL(tst_obj.vlidort_version_number, "2.8");
  BOOST_CHECK_EQUAL(tst_obj.vlidort_inunit, 21);
  BOOST_CHECK_EQUAL(tst_obj.vlidort_scenunit, 22);
  BOOST_CHECK_EQUAL(tst_obj.vlidort_funit, 23);
  BOOST_CHECK_EQUAL(tst_obj.vlidort_resunit, 24);
  BOOST_CHECK_EQUAL(tst_obj.vlidort_errunit, 25);
  BOOST_CHECK_EQUAL(tst_obj.maxstreams, 10);
  BOOST_CHECK_EQUAL(tst_obj.maxlayers, VLIDORT_MAXLAYER);
  BOOST_CHECK_EQUAL(tst_obj.maxfinelayers, 4);
  BOOST_CHECK_EQUAL(tst_obj.maxmoments_input, 200);
  BOOST_CHECK_EQUAL(tst_obj.max_thermal_coeffs, 2);
  BOOST_CHECK_EQUAL(tst_obj.max_szangles, 8);
  BOOST_CHECK_EQUAL(tst_obj.max_user_vzangles, 1);
  BOOST_CHECK_EQUAL(tst_obj.max_user_relazms, 1);
  BOOST_CHECK_EQUAL(tst_obj.max_user_obsgeoms, 8);
  BOOST_CHECK_EQUAL(tst_obj.max_user_levels, 2);
  BOOST_CHECK_EQUAL(tst_obj.max_partlayers, 2);
  BOOST_CHECK_EQUAL(tst_obj.max_taylor_terms, 7);
  BOOST_CHECK_EQUAL(tst_obj.maxstokes, 4);
  BOOST_CHECK_EQUAL(tst_obj.max_directions, 2);
  BOOST_CHECK_EQUAL(tst_obj.max_brdf_kernels, 4);
  BOOST_CHECK_EQUAL(tst_obj.max_brdf_parameters, 4);
  BOOST_CHECK_EQUAL(tst_obj.maxstreams_brdf, 100);
  BOOST_CHECK_EQUAL(tst_obj.max_msrs_muquad, 50);
  BOOST_CHECK_EQUAL(tst_obj.max_msrs_phiquad, 100);
  BOOST_CHECK_EQUAL(tst_obj.maxstreams_scaling, 24);
  BOOST_CHECK_EQUAL(tst_obj.max_atmoswfs, VLIDORT_MAXATMOSWFS);
  BOOST_CHECK_EQUAL(tst_obj.max_surfacewfs, 1);
  BOOST_CHECK_EQUAL(tst_obj.max_sleavewfs, 1);
  BOOST_CHECK_EQUAL(tst_obj.max_messages, 25);
  BOOST_CHECK_EQUAL(tst_obj.maxbeams, tst_obj.max_szangles);
  BOOST_CHECK_EQUAL(tst_obj.max_user_streams, tst_obj.max_user_vzangles);
  BOOST_CHECK_EQUAL(tst_obj.max_geometries, tst_obj.max_user_vzangles*tst_obj.max_user_relazms*tst_obj.max_szangles);
  BOOST_CHECK_EQUAL(tst_obj.max_allstrms, tst_obj.max_user_streams + tst_obj.maxstreams);
  BOOST_CHECK_EQUAL(tst_obj.max_allstrms_p1, tst_obj.max_allstrms + tst_obj.maxbeams);
  BOOST_CHECK_EQUAL(tst_obj.maxmoments, 2*tst_obj.maxstreams);
  BOOST_CHECK_EQUAL(tst_obj.maxfourier, 2*tst_obj.maxstreams - 1);
  BOOST_CHECK_EQUAL(tst_obj.maxstokes_sq, 16);
  BOOST_CHECK_EQUAL(tst_obj.maxsthalf_brdf, tst_obj.maxstreams_brdf / 2);
  BOOST_CHECK_EQUAL(tst_obj.maxstreams_2, 2*tst_obj.maxstreams);
  BOOST_CHECK_EQUAL(tst_obj.maxstreams_21, 2*tst_obj.maxstreams - 1);
  BOOST_CHECK_EQUAL(tst_obj.maxstreams_p1, tst_obj.maxstreams + 1);
  BOOST_CHECK_EQUAL(tst_obj.maxstreams_p2, tst_obj.maxstreams + 2);
  BOOST_CHECK_EQUAL(tst_obj.maxstrmstks, tst_obj.maxstreams * tst_obj.maxstokes);
  BOOST_CHECK_EQUAL(tst_obj.maxstrmstks_2, 2*tst_obj.maxstrmstks);
  BOOST_CHECK_EQUAL(tst_obj.maxstrmstks_21, 2*tst_obj.maxstrmstks - 1);
  BOOST_CHECK_EQUAL(tst_obj.maxstrmstks_p1, tst_obj.maxstrmstks + 1);
  BOOST_CHECK_EQUAL(tst_obj.maxstrmstks_p2, tst_obj.maxstrmstks + 2);
  BOOST_CHECK_EQUAL(tst_obj.maxstrmstks_p4, tst_obj.maxstrmstks + 4);
  BOOST_CHECK_EQUAL(tst_obj.max_ustrmstks, tst_obj.max_user_streams * tst_obj.maxstokes);
  BOOST_CHECK_EQUAL(tst_obj.maxevalues, tst_obj.maxstrmstks);
  BOOST_CHECK_EQUAL(tst_obj.maxtotal, tst_obj.maxlayers*tst_obj.maxstrmstks_2);
  BOOST_CHECK_EQUAL(tst_obj.maxbandtotal, 9*tst_obj.maxstrmstks - 2);
  BOOST_CHECK_EQUAL(tst_obj.max_psols, 2);
  BOOST_CHECK_EQUAL(tst_obj.max_scatpsols, tst_obj.max_psols);
  BOOST_CHECK_EQUAL(tst_obj.one, 1.0e0);
  BOOST_CHECK_EQUAL(tst_obj.zero, 0.0e0);
  BOOST_CHECK_EQUAL(tst_obj.onep5, 1.5e0);
  BOOST_CHECK_EQUAL(tst_obj.two, 2.0e0);
  BOOST_CHECK_EQUAL(tst_obj.three, 3.0e0);
  BOOST_CHECK_EQUAL(tst_obj.four, 4.0e0);
  BOOST_CHECK_EQUAL(tst_obj.quarter, 0.25e0);
  BOOST_CHECK_EQUAL(tst_obj.half, 0.5e0);
  BOOST_CHECK_EQUAL(tst_obj.minus_one, -tst_obj.one);
  BOOST_CHECK_EQUAL(tst_obj.minus_two, -tst_obj.two);
  BOOST_CHECK_EQUAL(tst_obj.deg_to_rad, 1.7453292519943d-02);
  BOOST_CHECK_EQUAL(tst_obj.pie, 180.0d0*tst_obj.deg_to_rad);
  BOOST_CHECK_EQUAL(tst_obj.pi2, 2.0d0 * tst_obj.pie);
  BOOST_CHECK_EQUAL(tst_obj.pi4, 4.0d0 * tst_obj.pie);
  BOOST_CHECK_EQUAL(tst_obj.pio2, tst_obj.half * tst_obj.pie);
  BOOST_CHECK_EQUAL(tst_obj.pio4, tst_obj.quarter * tst_obj.pie);
  BOOST_CHECK_EQUAL(tst_obj.eps3, 0.001e0);
  BOOST_CHECK_EQUAL(tst_obj.eps4, 0.0001e0);
  BOOST_CHECK_EQUAL(tst_obj.eps5, 0.00001e0);
  BOOST_CHECK_EQUAL(tst_obj.smallnum, 1.0d-9);
  BOOST_CHECK_EQUAL(tst_obj.bigexp, 32.0e0);
  BOOST_CHECK_EQUAL(tst_obj.taylor_small, 1.0d-04);
  BOOST_CHECK_EQUAL(tst_obj.taylor_large, 1.0d+04);
  BOOST_CHECK_EQUAL(tst_obj.hopital_tolerance, tst_obj.eps3);
  BOOST_CHECK_EQUAL(tst_obj.omega_smallnum, 1.0d-15);
  BOOST_CHECK_EQUAL(tst_obj.max_tau_spath, 32.0e0);
  BOOST_CHECK_EQUAL(tst_obj.max_tau_upath, 32.0e0);
  BOOST_CHECK_EQUAL(tst_obj.max_tau_qpath, 32.0e0);
  BOOST_CHECK_EQUAL(tst_obj.vlidort_serious, 4);
  BOOST_CHECK_EQUAL(tst_obj.vlidort_warning, 3);
  BOOST_CHECK_EQUAL(tst_obj.vlidort_info, 2);
  BOOST_CHECK_EQUAL(tst_obj.vlidort_debug, 1);
  BOOST_CHECK_EQUAL(tst_obj.vlidort_success, 0);
  BOOST_CHECK_EQUAL(tst_obj.upidx, 1);
  BOOST_CHECK_EQUAL(tst_obj.dnidx, 2);
  BOOST_CHECK_EQUAL(tst_obj.lambertian_idx, 1);
  BOOST_CHECK_EQUAL(tst_obj.rossthin_idx, 2);
  BOOST_CHECK_EQUAL(tst_obj.rossthick_idx, 3);
  BOOST_CHECK_EQUAL(tst_obj.lisparse_idx, 4);
  BOOST_CHECK_EQUAL(tst_obj.lidense_idx, 5);
  BOOST_CHECK_EQUAL(tst_obj.hapke_idx, 6);
  BOOST_CHECK_EQUAL(tst_obj.roujean_idx, 7);
  BOOST_CHECK_EQUAL(tst_obj.rahman_idx, 8);
  BOOST_CHECK_EQUAL(tst_obj.coxmunk_idx, 9);
  BOOST_CHECK_EQUAL(tst_obj.gisscoxmunk_idx, 10);
  BOOST_CHECK_EQUAL(tst_obj.gisscoxmunk_cri_idx, 11);
  BOOST_CHECK_EQUAL(tst_obj.bpdfsoil_idx, 12);
  BOOST_CHECK_EQUAL(tst_obj.bpdfvegn_idx, 13);
  BOOST_CHECK_EQUAL(tst_obj.bpdfndvi_idx, 14);
  BOOST_CHECK_EQUAL(tst_obj.newcmglint_idx, 15);
  BOOST_CHECK_EQUAL(tst_obj.newgcmglint_idx, 16);
  BOOST_CHECK_EQUAL(tst_obj.rtkhotspot_idx, 17);
  BOOST_CHECK_EQUAL(tst_obj.modfresnel_idx, 18);
  BOOST_CHECK_EQUAL(tst_obj.snowbrdf_idx, 19);
  BOOST_CHECK_EQUAL(tst_obj.maxbrdf_idx, tst_obj.snowbrdf_idx);
  
}

BOOST_AUTO_TEST_CASE(vbrdf_linsup_inputs)
{
  // Used for checking dimensions
  VLidort_Pars lid_pars = VLidort_Pars::instance();

  // Test constructor
  VBrdf_Linsup_Inputs tst_obj = VBrdf_Linsup_Inputs();

  // Test memory space sizes
  tst_obj.check_byte_sizes(); // throws exception on mismatch

  // Test variable shapes
  BOOST_CHECK_EQUAL(tst_obj.bs_do_kernel_factor_wfs().extent(0), lid_pars.max_brdf_kernels);
  BOOST_CHECK_EQUAL(tst_obj.bs_do_kernel_params_wfs().extent(0), lid_pars.max_brdf_kernels);
  BOOST_CHECK_EQUAL(tst_obj.bs_do_kernel_params_wfs().extent(1), lid_pars.max_brdf_parameters);
  BOOST_CHECK_EQUAL(tst_obj.bs_do_kparams_derivs().extent(0), lid_pars.max_brdf_kernels);
  

  // Test initialization
  BOOST_CHECK_EQUAL(tst_obj.bs_n_surface_wfs(), 0);
  BOOST_CHECK_EQUAL(tst_obj.bs_n_kernel_factor_wfs(), 0);
  BOOST_CHECK_EQUAL(tst_obj.bs_n_kernel_params_wfs(), 0);
  blitz::Array<bool, 1> bs_do_kernel_factor_wfs_exp(tst_obj.bs_do_kernel_factor_wfs().shape());
  bs_do_kernel_factor_wfs_exp = false;
  BOOST_CHECK_MATRIX_CLOSE_TOL(tst_obj.bs_do_kernel_factor_wfs(), bs_do_kernel_factor_wfs_exp, 1e-10);
  blitz::Array<bool, 2> bs_do_kernel_params_wfs_exp(tst_obj.bs_do_kernel_params_wfs().shape());
  bs_do_kernel_params_wfs_exp = false;
  BOOST_CHECK_MATRIX_CLOSE_TOL(tst_obj.bs_do_kernel_params_wfs(), bs_do_kernel_params_wfs_exp, 1e-10);
  blitz::Array<bool, 1> bs_do_kparams_derivs_exp(tst_obj.bs_do_kparams_derivs().shape());
  bs_do_kparams_derivs_exp = false;
  BOOST_CHECK_MATRIX_CLOSE_TOL(tst_obj.bs_do_kparams_derivs(), bs_do_kparams_derivs_exp, 1e-10);
  BOOST_CHECK_EQUAL(tst_obj.bs_do_bsavalue_wf(), false);
  BOOST_CHECK_EQUAL(tst_obj.bs_do_wsavalue_wf(), false);
  BOOST_CHECK_EQUAL(tst_obj.bs_do_windspeed_wf(), false);
  
}

BOOST_AUTO_TEST_CASE(vbrdf_linsup_outputs)
{
  // Used for checking dimensions
  VLidort_Pars lid_pars = VLidort_Pars::instance();

  // Test constructor
  VBrdf_Linsup_Outputs tst_obj = VBrdf_Linsup_Outputs();

  // Test memory space sizes
  tst_obj.check_byte_sizes(); // throws exception on mismatch

  // Test variable shapes
  BOOST_CHECK_EQUAL(tst_obj.bs_ls_dbounce_brdfunc().extent(0), lid_pars.max_surfacewfs);
  BOOST_CHECK_EQUAL(tst_obj.bs_ls_dbounce_brdfunc().extent(1), lid_pars.maxstokes_sq);
  BOOST_CHECK_EQUAL(tst_obj.bs_ls_dbounce_brdfunc().extent(2), lid_pars.max_user_streams);
  BOOST_CHECK_EQUAL(tst_obj.bs_ls_dbounce_brdfunc().extent(3), lid_pars.max_user_relazms);
  BOOST_CHECK_EQUAL(tst_obj.bs_ls_dbounce_brdfunc().extent(4), lid_pars.maxbeams);
  BOOST_CHECK_EQUAL(tst_obj.bs_ls_brdf_f_0().extent(0), lid_pars.max_surfacewfs);
  BOOST_CHECK_EQUAL(tst_obj.bs_ls_brdf_f_0().extent(1), lid_pars.maxmoments+1);
  BOOST_CHECK_EQUAL(tst_obj.bs_ls_brdf_f_0().extent(2), lid_pars.maxstokes_sq);
  BOOST_CHECK_EQUAL(tst_obj.bs_ls_brdf_f_0().extent(3), lid_pars.maxstreams);
  BOOST_CHECK_EQUAL(tst_obj.bs_ls_brdf_f_0().extent(4), lid_pars.maxbeams);
  BOOST_CHECK_EQUAL(tst_obj.bs_ls_brdf_f().extent(0), lid_pars.max_surfacewfs);
  BOOST_CHECK_EQUAL(tst_obj.bs_ls_brdf_f().extent(1), lid_pars.maxmoments+1);
  BOOST_CHECK_EQUAL(tst_obj.bs_ls_brdf_f().extent(2), lid_pars.maxstokes_sq);
  BOOST_CHECK_EQUAL(tst_obj.bs_ls_brdf_f().extent(3), lid_pars.maxstreams);
  BOOST_CHECK_EQUAL(tst_obj.bs_ls_brdf_f().extent(4), lid_pars.maxstreams);
  BOOST_CHECK_EQUAL(tst_obj.bs_ls_user_brdf_f_0().extent(0), lid_pars.max_surfacewfs);
  BOOST_CHECK_EQUAL(tst_obj.bs_ls_user_brdf_f_0().extent(1), lid_pars.maxmoments+1);
  BOOST_CHECK_EQUAL(tst_obj.bs_ls_user_brdf_f_0().extent(2), lid_pars.maxstokes_sq);
  BOOST_CHECK_EQUAL(tst_obj.bs_ls_user_brdf_f_0().extent(3), lid_pars.max_user_streams);
  BOOST_CHECK_EQUAL(tst_obj.bs_ls_user_brdf_f_0().extent(4), lid_pars.maxbeams);
  BOOST_CHECK_EQUAL(tst_obj.bs_ls_user_brdf_f().extent(0), lid_pars.max_surfacewfs);
  BOOST_CHECK_EQUAL(tst_obj.bs_ls_user_brdf_f().extent(1), lid_pars.maxmoments+1);
  BOOST_CHECK_EQUAL(tst_obj.bs_ls_user_brdf_f().extent(2), lid_pars.maxstokes_sq);
  BOOST_CHECK_EQUAL(tst_obj.bs_ls_user_brdf_f().extent(3), lid_pars.max_user_streams);
  BOOST_CHECK_EQUAL(tst_obj.bs_ls_user_brdf_f().extent(4), lid_pars.maxstreams);
  BOOST_CHECK_EQUAL(tst_obj.bs_ls_emissivity().extent(0), lid_pars.max_surfacewfs);
  BOOST_CHECK_EQUAL(tst_obj.bs_ls_emissivity().extent(1), lid_pars.maxstokes);
  BOOST_CHECK_EQUAL(tst_obj.bs_ls_emissivity().extent(2), lid_pars.maxstreams);
  BOOST_CHECK_EQUAL(tst_obj.bs_ls_user_emissivity().extent(0), lid_pars.max_surfacewfs);
  BOOST_CHECK_EQUAL(tst_obj.bs_ls_user_emissivity().extent(1), lid_pars.maxstokes);
  BOOST_CHECK_EQUAL(tst_obj.bs_ls_user_emissivity().extent(2), lid_pars.max_user_streams);
  

  // Test initialization
  blitz::Array<double, 5> bs_ls_dbounce_brdfunc_exp(tst_obj.bs_ls_dbounce_brdfunc().shape());
  bs_ls_dbounce_brdfunc_exp = 0;
  BOOST_CHECK_MATRIX_CLOSE_TOL(tst_obj.bs_ls_dbounce_brdfunc(), bs_ls_dbounce_brdfunc_exp, 1e-10);
  blitz::Array<double, 5> bs_ls_brdf_f_0_exp(tst_obj.bs_ls_brdf_f_0().shape());
  bs_ls_brdf_f_0_exp = 0;
  BOOST_CHECK_MATRIX_CLOSE_TOL(tst_obj.bs_ls_brdf_f_0(), bs_ls_brdf_f_0_exp, 1e-10);
  blitz::Array<double, 5> bs_ls_brdf_f_exp(tst_obj.bs_ls_brdf_f().shape());
  bs_ls_brdf_f_exp = 0;
  BOOST_CHECK_MATRIX_CLOSE_TOL(tst_obj.bs_ls_brdf_f(), bs_ls_brdf_f_exp, 1e-10);
  blitz::Array<double, 5> bs_ls_user_brdf_f_0_exp(tst_obj.bs_ls_user_brdf_f_0().shape());
  bs_ls_user_brdf_f_0_exp = 0;
  BOOST_CHECK_MATRIX_CLOSE_TOL(tst_obj.bs_ls_user_brdf_f_0(), bs_ls_user_brdf_f_0_exp, 1e-10);
  blitz::Array<double, 5> bs_ls_user_brdf_f_exp(tst_obj.bs_ls_user_brdf_f().shape());
  bs_ls_user_brdf_f_exp = 0;
  BOOST_CHECK_MATRIX_CLOSE_TOL(tst_obj.bs_ls_user_brdf_f(), bs_ls_user_brdf_f_exp, 1e-10);
  blitz::Array<double, 3> bs_ls_emissivity_exp(tst_obj.bs_ls_emissivity().shape());
  bs_ls_emissivity_exp = 0;
  BOOST_CHECK_MATRIX_CLOSE_TOL(tst_obj.bs_ls_emissivity(), bs_ls_emissivity_exp, 1e-10);
  blitz::Array<double, 3> bs_ls_user_emissivity_exp(tst_obj.bs_ls_user_emissivity().shape());
  bs_ls_user_emissivity_exp = 0;
  BOOST_CHECK_MATRIX_CLOSE_TOL(tst_obj.bs_ls_user_emissivity(), bs_ls_user_emissivity_exp, 1e-10);
  
}

BOOST_AUTO_TEST_CASE(vbrdf_sup_inputs)
{
  // Used for checking dimensions
  VLidort_Pars lid_pars = VLidort_Pars::instance();

  // Test constructor
  VBrdf_Sup_Inputs tst_obj = VBrdf_Sup_Inputs();

  // Test memory space sizes
  tst_obj.check_byte_sizes(); // throws exception on mismatch

  // Test variable shapes
  BOOST_CHECK_EQUAL(tst_obj.bs_beam_szas().extent(0), lid_pars.maxbeams);
  BOOST_CHECK_EQUAL(tst_obj.bs_user_angles_input().extent(0), lid_pars.max_user_streams);
  BOOST_CHECK_EQUAL(tst_obj.bs_user_relazms().extent(0), lid_pars.max_user_relazms);
  BOOST_CHECK_EQUAL(tst_obj.bs_user_obsgeoms().extent(0), lid_pars.max_user_obsgeoms);
  BOOST_CHECK_EQUAL(tst_obj.bs_user_obsgeoms().extent(1), lid_pars.3);
  BOOST_CHECK_EQUAL(tst_obj.bs_user_doublets().extent(0), lid_pars.max_user_streams);
  BOOST_CHECK_EQUAL(tst_obj.bs_user_doublets().extent(1), lid_pars.2);
  BOOST_CHECK_EQUAL(tst_obj.bs_which_brdf().extent(0), lid_pars.max_brdf_kernels);
  BOOST_CHECK_EQUAL(tst_obj.bs_lambertian_kernel_flag().extent(0), lid_pars.max_brdf_kernels);
  BOOST_CHECK_EQUAL(tst_obj.bs_brdf_factors().extent(0), lid_pars.max_brdf_kernels);
  BOOST_CHECK_EQUAL(tst_obj.bs_n_brdf_parameters().extent(0), lid_pars.max_brdf_kernels);
  BOOST_CHECK_EQUAL(tst_obj.bs_brdf_parameters().extent(0), lid_pars.max_brdf_kernels);
  BOOST_CHECK_EQUAL(tst_obj.bs_brdf_parameters().extent(1), lid_pars.max_brdf_parameters);
  BOOST_CHECK_EQUAL(tst_obj.bs_winddir().extent(0), lid_pars.maxbeams);
  

  // Test initialization
  BOOST_CHECK_EQUAL(tst_obj.bs_do_brdf_surface(), false);
  BOOST_CHECK_EQUAL(tst_obj.bs_do_user_streams(), false);
  BOOST_CHECK_EQUAL(tst_obj.bs_do_solar_sources(), false);
  BOOST_CHECK_EQUAL(tst_obj.bs_do_user_obsgeoms(), false);
  BOOST_CHECK_EQUAL(tst_obj.bs_do_doublet_geometry(), false);
  BOOST_CHECK_EQUAL(tst_obj.bs_do_surface_emission(), false);
  BOOST_CHECK_EQUAL(tst_obj.bs_nstokes(), 0);
  BOOST_CHECK_EQUAL(tst_obj.bs_nstreams(), 0);
  BOOST_CHECK_EQUAL(tst_obj.bs_nbeams(), 0);
  blitz::Array<double, 1> bs_beam_szas_exp(tst_obj.bs_beam_szas().shape());
  bs_beam_szas_exp = 0;
  BOOST_CHECK_MATRIX_CLOSE_TOL(tst_obj.bs_beam_szas(), bs_beam_szas_exp, 1e-10);
  BOOST_CHECK_EQUAL(tst_obj.bs_n_user_streams(), 0);
  blitz::Array<double, 1> bs_user_angles_input_exp(tst_obj.bs_user_angles_input().shape());
  bs_user_angles_input_exp = 0;
  BOOST_CHECK_MATRIX_CLOSE_TOL(tst_obj.bs_user_angles_input(), bs_user_angles_input_exp, 1e-10);
  BOOST_CHECK_EQUAL(tst_obj.bs_n_user_relazms(), 0);
  blitz::Array<double, 1> bs_user_relazms_exp(tst_obj.bs_user_relazms().shape());
  bs_user_relazms_exp = 0;
  BOOST_CHECK_MATRIX_CLOSE_TOL(tst_obj.bs_user_relazms(), bs_user_relazms_exp, 1e-10);
  BOOST_CHECK_EQUAL(tst_obj.bs_n_user_obsgeoms(), 0);
  blitz::Array<double, 2> bs_user_obsgeoms_exp(tst_obj.bs_user_obsgeoms().shape());
  bs_user_obsgeoms_exp = 0;
  BOOST_CHECK_MATRIX_CLOSE_TOL(tst_obj.bs_user_obsgeoms(), bs_user_obsgeoms_exp, 1e-10);
  BOOST_CHECK_EQUAL(tst_obj.bs_n_user_doublets(), 0);
  blitz::Array<double, 2> bs_user_doublets_exp(tst_obj.bs_user_doublets().shape());
  bs_user_doublets_exp = 0;
  BOOST_CHECK_MATRIX_CLOSE_TOL(tst_obj.bs_user_doublets(), bs_user_doublets_exp, 1e-10);
  BOOST_CHECK_EQUAL(tst_obj.bs_nstreams_brdf(), 0);
  BOOST_CHECK_EQUAL(tst_obj.bs_n_brdf_kernels(), 0);
  blitz::Array<int, 1> bs_which_brdf_exp(tst_obj.bs_which_brdf().shape());
  bs_which_brdf_exp = 0;
  BOOST_CHECK_MATRIX_CLOSE_TOL(tst_obj.bs_which_brdf(), bs_which_brdf_exp, 1e-10);
  blitz::Array<bool, 1> bs_lambertian_kernel_flag_exp(tst_obj.bs_lambertian_kernel_flag().shape());
  bs_lambertian_kernel_flag_exp = false;
  BOOST_CHECK_MATRIX_CLOSE_TOL(tst_obj.bs_lambertian_kernel_flag(), bs_lambertian_kernel_flag_exp, 1e-10);
  blitz::Array<double, 1> bs_brdf_factors_exp(tst_obj.bs_brdf_factors().shape());
  bs_brdf_factors_exp = 0;
  BOOST_CHECK_MATRIX_CLOSE_TOL(tst_obj.bs_brdf_factors(), bs_brdf_factors_exp, 1e-10);
  blitz::Array<int, 1> bs_n_brdf_parameters_exp(tst_obj.bs_n_brdf_parameters().shape());
  bs_n_brdf_parameters_exp = 0;
  BOOST_CHECK_MATRIX_CLOSE_TOL(tst_obj.bs_n_brdf_parameters(), bs_n_brdf_parameters_exp, 1e-10);
  blitz::Array<double, 2> bs_brdf_parameters_exp(tst_obj.bs_brdf_parameters().shape());
  bs_brdf_parameters_exp = 0;
  BOOST_CHECK_MATRIX_CLOSE_TOL(tst_obj.bs_brdf_parameters(), bs_brdf_parameters_exp, 1e-10);
  BOOST_CHECK_EQUAL(tst_obj.bs_do_directbounce_only(), false);
  BOOST_CHECK_EQUAL(tst_obj.bs_do_shadow_effect(), false);
  BOOST_CHECK_EQUAL(tst_obj.bs_do_wsabsa_output(), false);
  BOOST_CHECK_EQUAL(tst_obj.bs_do_wsa_scaling(), false);
  BOOST_CHECK_EQUAL(tst_obj.bs_do_bsa_scaling(), false);
  BOOST_CHECK_CLOSE(tst_obj.bs_wsa_value(), 0, 1e-10);
  BOOST_CHECK_CLOSE(tst_obj.bs_bsa_value(), 0, 1e-10);
  BOOST_CHECK_EQUAL(tst_obj.bs_do_newcmglint(), false);
  BOOST_CHECK_EQUAL(tst_obj.bs_do_newgcmglint(), false);
  BOOST_CHECK_CLOSE(tst_obj.bs_salinity(), 0, 1e-10);
  BOOST_CHECK_CLOSE(tst_obj.bs_wavelength(), 0, 1e-10);
  BOOST_CHECK_CLOSE(tst_obj.bs_windspeed(), 0, 1e-10);
  blitz::Array<double, 1> bs_winddir_exp(tst_obj.bs_winddir().shape());
  bs_winddir_exp = 0;
  BOOST_CHECK_MATRIX_CLOSE_TOL(tst_obj.bs_winddir(), bs_winddir_exp, 1e-10);
  BOOST_CHECK_EQUAL(tst_obj.bs_do_glintshadow(), false);
  BOOST_CHECK_EQUAL(tst_obj.bs_do_foamoption(), false);
  BOOST_CHECK_EQUAL(tst_obj.bs_do_facetisotropy(), false);
  BOOST_CHECK_EQUAL(tst_obj.bs_do_glitter_msrcorr(), false);
  BOOST_CHECK_EQUAL(tst_obj.bs_do_glitter_msrcorr_dbonly(), false);
  BOOST_CHECK_EQUAL(tst_obj.bs_glitter_msrcorr_order(), 0);
  BOOST_CHECK_EQUAL(tst_obj.bs_glitter_msrcorr_nmuquad(), 0);
  BOOST_CHECK_EQUAL(tst_obj.bs_glitter_msrcorr_nphiquad(), 0);
  
}

BOOST_AUTO_TEST_CASE(vbrdf_sup_outputs)
{
  // Used for checking dimensions
  VLidort_Pars lid_pars = VLidort_Pars::instance();

  // Test constructor
  VBrdf_Sup_Outputs tst_obj = VBrdf_Sup_Outputs();

  // Test memory space sizes
  tst_obj.check_byte_sizes(); // throws exception on mismatch

  // Test variable shapes
  BOOST_CHECK_EQUAL(tst_obj.bs_dbounce_brdfunc().extent(0), lid_pars.maxstokes_sq);
  BOOST_CHECK_EQUAL(tst_obj.bs_dbounce_brdfunc().extent(1), lid_pars.max_user_streams);
  BOOST_CHECK_EQUAL(tst_obj.bs_dbounce_brdfunc().extent(2), lid_pars.max_user_relazms);
  BOOST_CHECK_EQUAL(tst_obj.bs_dbounce_brdfunc().extent(3), lid_pars.maxbeams);
  BOOST_CHECK_EQUAL(tst_obj.bs_brdf_f_0().extent(0), lid_pars.maxmoments+1);
  BOOST_CHECK_EQUAL(tst_obj.bs_brdf_f_0().extent(1), lid_pars.maxstokes_sq);
  BOOST_CHECK_EQUAL(tst_obj.bs_brdf_f_0().extent(2), lid_pars.maxstreams);
  BOOST_CHECK_EQUAL(tst_obj.bs_brdf_f_0().extent(3), lid_pars.maxbeams);
  BOOST_CHECK_EQUAL(tst_obj.bs_brdf_f().extent(0), lid_pars.maxmoments+1);
  BOOST_CHECK_EQUAL(tst_obj.bs_brdf_f().extent(1), lid_pars.maxstokes_sq);
  BOOST_CHECK_EQUAL(tst_obj.bs_brdf_f().extent(2), lid_pars.maxstreams);
  BOOST_CHECK_EQUAL(tst_obj.bs_brdf_f().extent(3), lid_pars.maxstreams);
  BOOST_CHECK_EQUAL(tst_obj.bs_user_brdf_f_0().extent(0), lid_pars.maxmoments+1);
  BOOST_CHECK_EQUAL(tst_obj.bs_user_brdf_f_0().extent(1), lid_pars.maxstokes_sq);
  BOOST_CHECK_EQUAL(tst_obj.bs_user_brdf_f_0().extent(2), lid_pars.max_user_streams);
  BOOST_CHECK_EQUAL(tst_obj.bs_user_brdf_f_0().extent(3), lid_pars.maxbeams);
  BOOST_CHECK_EQUAL(tst_obj.bs_user_brdf_f().extent(0), lid_pars.maxmoments+1);
  BOOST_CHECK_EQUAL(tst_obj.bs_user_brdf_f().extent(1), lid_pars.maxstokes_sq);
  BOOST_CHECK_EQUAL(tst_obj.bs_user_brdf_f().extent(2), lid_pars.max_user_streams);
  BOOST_CHECK_EQUAL(tst_obj.bs_user_brdf_f().extent(3), lid_pars.maxstreams);
  BOOST_CHECK_EQUAL(tst_obj.bs_emissivity().extent(0), lid_pars.maxstokes);
  BOOST_CHECK_EQUAL(tst_obj.bs_emissivity().extent(1), lid_pars.maxstreams);
  BOOST_CHECK_EQUAL(tst_obj.bs_user_emissivity().extent(0), lid_pars.maxstokes);
  BOOST_CHECK_EQUAL(tst_obj.bs_user_emissivity().extent(1), lid_pars.max_user_streams);
  BOOST_CHECK_EQUAL(tst_obj.bs_wsa_kernels().extent(0), lid_pars.max_brdf_kernels);
  BOOST_CHECK_EQUAL(tst_obj.bs_bsa_kernels().extent(0), lid_pars.max_brdf_kernels);
  

  // Test initialization
  blitz::Array<double, 4> bs_dbounce_brdfunc_exp(tst_obj.bs_dbounce_brdfunc().shape());
  bs_dbounce_brdfunc_exp = 0;
  BOOST_CHECK_MATRIX_CLOSE_TOL(tst_obj.bs_dbounce_brdfunc(), bs_dbounce_brdfunc_exp, 1e-10);
  blitz::Array<double, 4> bs_brdf_f_0_exp(tst_obj.bs_brdf_f_0().shape());
  bs_brdf_f_0_exp = 0;
  BOOST_CHECK_MATRIX_CLOSE_TOL(tst_obj.bs_brdf_f_0(), bs_brdf_f_0_exp, 1e-10);
  blitz::Array<double, 4> bs_brdf_f_exp(tst_obj.bs_brdf_f().shape());
  bs_brdf_f_exp = 0;
  BOOST_CHECK_MATRIX_CLOSE_TOL(tst_obj.bs_brdf_f(), bs_brdf_f_exp, 1e-10);
  blitz::Array<double, 4> bs_user_brdf_f_0_exp(tst_obj.bs_user_brdf_f_0().shape());
  bs_user_brdf_f_0_exp = 0;
  BOOST_CHECK_MATRIX_CLOSE_TOL(tst_obj.bs_user_brdf_f_0(), bs_user_brdf_f_0_exp, 1e-10);
  blitz::Array<double, 4> bs_user_brdf_f_exp(tst_obj.bs_user_brdf_f().shape());
  bs_user_brdf_f_exp = 0;
  BOOST_CHECK_MATRIX_CLOSE_TOL(tst_obj.bs_user_brdf_f(), bs_user_brdf_f_exp, 1e-10);
  blitz::Array<double, 2> bs_emissivity_exp(tst_obj.bs_emissivity().shape());
  bs_emissivity_exp = 0;
  BOOST_CHECK_MATRIX_CLOSE_TOL(tst_obj.bs_emissivity(), bs_emissivity_exp, 1e-10);
  blitz::Array<double, 2> bs_user_emissivity_exp(tst_obj.bs_user_emissivity().shape());
  bs_user_emissivity_exp = 0;
  BOOST_CHECK_MATRIX_CLOSE_TOL(tst_obj.bs_user_emissivity(), bs_user_emissivity_exp, 1e-10);
  BOOST_CHECK_CLOSE(tst_obj.bs_wsa_calculated(), 0, 1e-10);
  blitz::Array<double, 1> bs_wsa_kernels_exp(tst_obj.bs_wsa_kernels().shape());
  bs_wsa_kernels_exp = 0;
  BOOST_CHECK_MATRIX_CLOSE_TOL(tst_obj.bs_wsa_kernels(), bs_wsa_kernels_exp, 1e-10);
  BOOST_CHECK_CLOSE(tst_obj.bs_bsa_calculated(), 0, 1e-10);
  blitz::Array<double, 1> bs_bsa_kernels_exp(tst_obj.bs_bsa_kernels().shape());
  bs_bsa_kernels_exp = 0;
  BOOST_CHECK_MATRIX_CLOSE_TOL(tst_obj.bs_bsa_kernels(), bs_bsa_kernels_exp, 1e-10);
  
}

BOOST_AUTO_TEST_CASE(vbrdf_input_exception_handling)
{
  // Used for checking dimensions
  VLidort_Pars lid_pars = VLidort_Pars::instance();

  // Test constructor
  VBrdf_Input_Exception_Handling tst_obj = VBrdf_Input_Exception_Handling();

  // Test memory space sizes
  tst_obj.check_byte_sizes(); // throws exception on mismatch

  // Test variable shapes
  

  // Test initialization
  BOOST_CHECK_EQUAL(tst_obj.bs_status_inputread(), 0);
  BOOST_CHECK_EQUAL(tst_obj.bs_ninputmessages(), 0);
  
}

BOOST_AUTO_TEST_CASE(vbrdf_output_exception_handling)
{
  // Used for checking dimensions
  VLidort_Pars lid_pars = VLidort_Pars::instance();

  // Test constructor
  VBrdf_Output_Exception_Handling tst_obj = VBrdf_Output_Exception_Handling();

  // Test memory space sizes
  tst_obj.check_byte_sizes(); // throws exception on mismatch

  // Test variable shapes
  

  // Test initialization
  BOOST_CHECK_EQUAL(tst_obj.bs_status_output(), 0);
  BOOST_CHECK_EQUAL(tst_obj.bs_noutputmessages(), 0);
  
}

BOOST_AUTO_TEST_CASE(vsleave_sup_inputs)
{
  // Used for checking dimensions
  VLidort_Pars lid_pars = VLidort_Pars::instance();

  // Test constructor
  Vsleave_Sup_Inputs tst_obj = Vsleave_Sup_Inputs();

  // Test memory space sizes
  tst_obj.check_byte_sizes(); // throws exception on mismatch

  // Test variable shapes
  BOOST_CHECK_EQUAL(tst_obj.sl_beam_szas().extent(0), lid_pars.maxbeams);
  BOOST_CHECK_EQUAL(tst_obj.sl_user_relazms().extent(0), lid_pars.max_user_relazms);
  BOOST_CHECK_EQUAL(tst_obj.sl_user_angles_input().extent(0), lid_pars.max_user_streams);
  BOOST_CHECK_EQUAL(tst_obj.sl_user_obsgeoms().extent(0), lid_pars.max_user_obsgeoms);
  BOOST_CHECK_EQUAL(tst_obj.sl_user_obsgeoms().extent(1), lid_pars.3);
  BOOST_CHECK_EQUAL(tst_obj.sl_user_doublets().extent(0), lid_pars.max_user_streams);
  BOOST_CHECK_EQUAL(tst_obj.sl_user_doublets().extent(1), lid_pars.2);
  BOOST_CHECK_EQUAL(tst_obj.sl_winddir().extent(0), lid_pars.maxbeams);
  BOOST_CHECK_EQUAL(tst_obj.sl_fl_epoch().extent(0), lid_pars.6);
  BOOST_CHECK_EQUAL(tst_obj.sl_fl_inputgaussians().extent(0), lid_pars.3);
  BOOST_CHECK_EQUAL(tst_obj.sl_fl_inputgaussians().extent(1), lid_pars.2);
  

  // Test initialization
  BOOST_CHECK_EQUAL(tst_obj.sl_do_sleaving(), false);
  BOOST_CHECK_EQUAL(tst_obj.sl_do_isotropic(), false);
  BOOST_CHECK_EQUAL(tst_obj.sl_do_roughsurface(), false);
  BOOST_CHECK_EQUAL(tst_obj.sl_do_exact(), false);
  BOOST_CHECK_EQUAL(tst_obj.sl_do_exactonly(), false);
  BOOST_CHECK_EQUAL(tst_obj.sl_do_fluorescence(), false);
  BOOST_CHECK_EQUAL(tst_obj.sl_do_solar_sources(), false);
  BOOST_CHECK_EQUAL(tst_obj.sl_do_user_streams(), false);
  BOOST_CHECK_EQUAL(tst_obj.sl_do_user_obsgeoms(), false);
  BOOST_CHECK_EQUAL(tst_obj.sl_do_doublet_geometry(), false);
  BOOST_CHECK_EQUAL(tst_obj.sl_nstokes(), 0);
  BOOST_CHECK_EQUAL(tst_obj.sl_nstreams(), 0);
  BOOST_CHECK_EQUAL(tst_obj.sl_nbeams(), 0);
  blitz::Array<double, 1> sl_beam_szas_exp(tst_obj.sl_beam_szas().shape());
  sl_beam_szas_exp = 0;
  BOOST_CHECK_MATRIX_CLOSE_TOL(tst_obj.sl_beam_szas(), sl_beam_szas_exp, 1e-10);
  BOOST_CHECK_EQUAL(tst_obj.sl_n_user_relazms(), 0);
  blitz::Array<double, 1> sl_user_relazms_exp(tst_obj.sl_user_relazms().shape());
  sl_user_relazms_exp = 0;
  BOOST_CHECK_MATRIX_CLOSE_TOL(tst_obj.sl_user_relazms(), sl_user_relazms_exp, 1e-10);
  BOOST_CHECK_EQUAL(tst_obj.sl_n_user_streams(), 0);
  blitz::Array<double, 1> sl_user_angles_input_exp(tst_obj.sl_user_angles_input().shape());
  sl_user_angles_input_exp = 0;
  BOOST_CHECK_MATRIX_CLOSE_TOL(tst_obj.sl_user_angles_input(), sl_user_angles_input_exp, 1e-10);
  BOOST_CHECK_EQUAL(tst_obj.sl_n_user_obsgeoms(), 0);
  blitz::Array<double, 2> sl_user_obsgeoms_exp(tst_obj.sl_user_obsgeoms().shape());
  sl_user_obsgeoms_exp = 0;
  BOOST_CHECK_MATRIX_CLOSE_TOL(tst_obj.sl_user_obsgeoms(), sl_user_obsgeoms_exp, 1e-10);
  BOOST_CHECK_EQUAL(tst_obj.sl_n_user_doublets(), 0);
  blitz::Array<double, 2> sl_user_doublets_exp(tst_obj.sl_user_doublets().shape());
  sl_user_doublets_exp = 0;
  BOOST_CHECK_MATRIX_CLOSE_TOL(tst_obj.sl_user_doublets(), sl_user_doublets_exp, 1e-10);
  BOOST_CHECK_CLOSE(tst_obj.sl_salinity(), 0, 1e-10);
  BOOST_CHECK_CLOSE(tst_obj.sl_chlorconc(), 0, 1e-10);
  BOOST_CHECK_CLOSE(tst_obj.sl_wavelength(), 0, 1e-10);
  BOOST_CHECK_EQUAL(tst_obj.sl_azimuthdep(), false);
  BOOST_CHECK_EQUAL(tst_obj.sl_do_fourier_output(), false);
  BOOST_CHECK_CLOSE(tst_obj.sl_windspeed(), 0, 1e-10);
  blitz::Array<double, 1> sl_winddir_exp(tst_obj.sl_winddir().shape());
  sl_winddir_exp = 0;
  BOOST_CHECK_MATRIX_CLOSE_TOL(tst_obj.sl_winddir(), sl_winddir_exp, 1e-10);
  BOOST_CHECK_EQUAL(tst_obj.sl_do_glintshadow(), false);
  BOOST_CHECK_EQUAL(tst_obj.sl_do_foamoption(), false);
  BOOST_CHECK_EQUAL(tst_obj.sl_do_facetisotropy(), false);
  BOOST_CHECK_CLOSE(tst_obj.sl_fl_wavelength(), 0, 1e-10);
  BOOST_CHECK_CLOSE(tst_obj.sl_fl_latitude(), 0, 1e-10);
  BOOST_CHECK_CLOSE(tst_obj.sl_fl_longitude(), 0, 1e-10);
  blitz::Array<int, 1> sl_fl_epoch_exp(tst_obj.sl_fl_epoch().shape());
  sl_fl_epoch_exp = 0;
  BOOST_CHECK_MATRIX_CLOSE_TOL(tst_obj.sl_fl_epoch(), sl_fl_epoch_exp, 1e-10);
  BOOST_CHECK_CLOSE(tst_obj.sl_fl_amplitude755(), 0, 1e-10);
  BOOST_CHECK_EQUAL(tst_obj.sl_fl_do_datagaussian(), false);
  blitz::Array<double, 2> sl_fl_inputgaussians_exp(tst_obj.sl_fl_inputgaussians().shape());
  sl_fl_inputgaussians_exp = 0;
  BOOST_CHECK_MATRIX_CLOSE_TOL(tst_obj.sl_fl_inputgaussians(), sl_fl_inputgaussians_exp, 1e-10);
  
}

BOOST_AUTO_TEST_CASE(vlidort_fixed_lincontrol)
{
  // Used for checking dimensions
  VLidort_Pars lid_pars = VLidort_Pars::instance();

  // Test constructor
  VLidort_Fixed_Lincontrol tst_obj = VLidort_Fixed_Lincontrol();

  // Test memory space sizes
  tst_obj.check_byte_sizes(); // throws exception on mismatch

  // Test variable shapes
  BOOST_CHECK_EQUAL(tst_obj.ts_layer_vary_flag().extent(0), lid_pars.maxlayers);
  BOOST_CHECK_EQUAL(tst_obj.ts_layer_vary_number().extent(0), lid_pars.maxlayers);
  

  // Test initialization
  blitz::Array<bool, 1> ts_layer_vary_flag_exp(tst_obj.ts_layer_vary_flag().shape());
  ts_layer_vary_flag_exp = false;
  BOOST_CHECK_MATRIX_CLOSE_TOL(tst_obj.ts_layer_vary_flag(), ts_layer_vary_flag_exp, 1e-10);
  blitz::Array<int, 1> ts_layer_vary_number_exp(tst_obj.ts_layer_vary_number().shape());
  ts_layer_vary_number_exp = 0;
  BOOST_CHECK_MATRIX_CLOSE_TOL(tst_obj.ts_layer_vary_number(), ts_layer_vary_number_exp, 1e-10);
  BOOST_CHECK_EQUAL(tst_obj.ts_n_totalcolumn_wfs(), 0);
  BOOST_CHECK_EQUAL(tst_obj.ts_n_totalprofile_wfs(), 0);
  BOOST_CHECK_EQUAL(tst_obj.ts_n_surface_wfs(), 0);
  BOOST_CHECK_EQUAL(tst_obj.ts_n_sleave_wfs(), 0);
  
}

BOOST_AUTO_TEST_CASE(vlidort_fixed_linoptical)
{
  // Used for checking dimensions
  VLidort_Pars lid_pars = VLidort_Pars::instance();

  // Test constructor
  VLidort_Fixed_Linoptical tst_obj = VLidort_Fixed_Linoptical();

  // Test memory space sizes
  tst_obj.check_byte_sizes(); // throws exception on mismatch

  // Test variable shapes
  BOOST_CHECK_EQUAL(tst_obj.ts_l_deltau_vert_input().extent(0), lid_pars.max_atmoswfs);
  BOOST_CHECK_EQUAL(tst_obj.ts_l_deltau_vert_input().extent(1), lid_pars.maxlayers);
  BOOST_CHECK_EQUAL(tst_obj.ts_l_omega_total_input().extent(0), lid_pars.max_atmoswfs);
  BOOST_CHECK_EQUAL(tst_obj.ts_l_omega_total_input().extent(1), lid_pars.maxlayers);
  BOOST_CHECK_EQUAL(tst_obj.ts_l_greekmat_total_input().extent(0), lid_pars.max_atmoswfs);
  BOOST_CHECK_EQUAL(tst_obj.ts_l_greekmat_total_input().extent(1), lid_pars.maxmoments_input+1);
  BOOST_CHECK_EQUAL(tst_obj.ts_l_greekmat_total_input().extent(2), lid_pars.maxlayers);
  BOOST_CHECK_EQUAL(tst_obj.ts_l_greekmat_total_input().extent(3), lid_pars.maxstokes_sq);
  BOOST_CHECK_EQUAL(tst_obj.ts_l_fmatrix_up().extent(0), lid_pars.max_atmoswfs);
  BOOST_CHECK_EQUAL(tst_obj.ts_l_fmatrix_up().extent(1), lid_pars.maxlayers);
  BOOST_CHECK_EQUAL(tst_obj.ts_l_fmatrix_up().extent(2), lid_pars.max_geometries);
  BOOST_CHECK_EQUAL(tst_obj.ts_l_fmatrix_up().extent(3), lid_pars.6);
  BOOST_CHECK_EQUAL(tst_obj.ts_l_fmatrix_dn().extent(0), lid_pars.max_atmoswfs);
  BOOST_CHECK_EQUAL(tst_obj.ts_l_fmatrix_dn().extent(1), lid_pars.maxlayers);
  BOOST_CHECK_EQUAL(tst_obj.ts_l_fmatrix_dn().extent(2), lid_pars.max_geometries);
  BOOST_CHECK_EQUAL(tst_obj.ts_l_fmatrix_dn().extent(3), lid_pars.6);
  

  // Test initialization
  blitz::Array<double, 2> ts_l_deltau_vert_input_exp(tst_obj.ts_l_deltau_vert_input().shape());
  ts_l_deltau_vert_input_exp = 0;
  BOOST_CHECK_MATRIX_CLOSE_TOL(tst_obj.ts_l_deltau_vert_input(), ts_l_deltau_vert_input_exp, 1e-10);
  blitz::Array<double, 2> ts_l_omega_total_input_exp(tst_obj.ts_l_omega_total_input().shape());
  ts_l_omega_total_input_exp = 0;
  BOOST_CHECK_MATRIX_CLOSE_TOL(tst_obj.ts_l_omega_total_input(), ts_l_omega_total_input_exp, 1e-10);
  blitz::Array<double, 4> ts_l_greekmat_total_input_exp(tst_obj.ts_l_greekmat_total_input().shape());
  ts_l_greekmat_total_input_exp = 0;
  BOOST_CHECK_MATRIX_CLOSE_TOL(tst_obj.ts_l_greekmat_total_input(), ts_l_greekmat_total_input_exp, 1e-10);
  blitz::Array<double, 4> ts_l_fmatrix_up_exp(tst_obj.ts_l_fmatrix_up().shape());
  ts_l_fmatrix_up_exp = 0;
  BOOST_CHECK_MATRIX_CLOSE_TOL(tst_obj.ts_l_fmatrix_up(), ts_l_fmatrix_up_exp, 1e-10);
  blitz::Array<double, 4> ts_l_fmatrix_dn_exp(tst_obj.ts_l_fmatrix_dn().shape());
  ts_l_fmatrix_dn_exp = 0;
  BOOST_CHECK_MATRIX_CLOSE_TOL(tst_obj.ts_l_fmatrix_dn(), ts_l_fmatrix_dn_exp, 1e-10);
  
}

BOOST_AUTO_TEST_CASE(vlidort_fixed_lininputs)
{
  // Used for checking dimensions
  VLidort_Pars lid_pars = VLidort_Pars::instance();

  // Test constructor
  VLidort_Fixed_Lininputs tst_obj = VLidort_Fixed_Lininputs();

  // Test memory space sizes
  tst_obj.check_byte_sizes(); // throws exception on mismatch

  // Test variable shapes
  

  // Test initialization
  
}

BOOST_AUTO_TEST_CASE(vlidort_modified_lincontrol)
{
  // Used for checking dimensions
  VLidort_Pars lid_pars = VLidort_Pars::instance();

  // Test constructor
  VLidort_Modified_Lincontrol tst_obj = VLidort_Modified_Lincontrol();

  // Test memory space sizes
  tst_obj.check_byte_sizes(); // throws exception on mismatch

  // Test variable shapes
  

  // Test initialization
  BOOST_CHECK_EQUAL(tst_obj.ts_do_column_linearization(), false);
  BOOST_CHECK_EQUAL(tst_obj.ts_do_profile_linearization(), false);
  BOOST_CHECK_EQUAL(tst_obj.ts_do_atmos_linearization(), false);
  BOOST_CHECK_EQUAL(tst_obj.ts_do_surface_linearization(), false);
  BOOST_CHECK_EQUAL(tst_obj.ts_do_linearization(), false);
  BOOST_CHECK_EQUAL(tst_obj.ts_do_simulation_only(), false);
  BOOST_CHECK_EQUAL(tst_obj.ts_do_atmos_lbbf(), false);
  BOOST_CHECK_EQUAL(tst_obj.ts_do_surface_lbbf(), false);
  BOOST_CHECK_EQUAL(tst_obj.ts_do_sleave_wfs(), false);
  
}

BOOST_AUTO_TEST_CASE(vlidort_modified_lininputs)
{
  // Used for checking dimensions
  VLidort_Pars lid_pars = VLidort_Pars::instance();

  // Test constructor
  VLidort_Modified_Lininputs tst_obj = VLidort_Modified_Lininputs();

  // Test memory space sizes
  tst_obj.check_byte_sizes(); // throws exception on mismatch

  // Test variable shapes
  

  // Test initialization
  
}

BOOST_AUTO_TEST_CASE(vlidort_lincol)
{
  // Used for checking dimensions
  VLidort_Pars lid_pars = VLidort_Pars::instance();

  // Test constructor
  VLidort_Lincol tst_obj = VLidort_Lincol();

  // Test memory space sizes
  tst_obj.check_byte_sizes(); // throws exception on mismatch

  // Test variable shapes
  BOOST_CHECK_EQUAL(tst_obj.ts_columnwf().extent(0), lid_pars.max_atmoswfs);
  BOOST_CHECK_EQUAL(tst_obj.ts_columnwf().extent(1), lid_pars.max_user_levels);
  BOOST_CHECK_EQUAL(tst_obj.ts_columnwf().extent(2), lid_pars.max_geometries);
  BOOST_CHECK_EQUAL(tst_obj.ts_columnwf().extent(3), lid_pars.maxstokes);
  BOOST_CHECK_EQUAL(tst_obj.ts_columnwf().extent(4), lid_pars.max_directions);
  BOOST_CHECK_EQUAL(tst_obj.ts_meanst_diffuse_colwf().extent(0), lid_pars.max_atmoswfs);
  BOOST_CHECK_EQUAL(tst_obj.ts_meanst_diffuse_colwf().extent(1), lid_pars.max_user_levels);
  BOOST_CHECK_EQUAL(tst_obj.ts_meanst_diffuse_colwf().extent(2), lid_pars.max_szangles);
  BOOST_CHECK_EQUAL(tst_obj.ts_meanst_diffuse_colwf().extent(3), lid_pars.maxstokes);
  BOOST_CHECK_EQUAL(tst_obj.ts_meanst_diffuse_colwf().extent(4), lid_pars.max_directions);
  BOOST_CHECK_EQUAL(tst_obj.ts_flux_diffuse_colwf().extent(0), lid_pars.max_atmoswfs);
  BOOST_CHECK_EQUAL(tst_obj.ts_flux_diffuse_colwf().extent(1), lid_pars.max_user_levels);
  BOOST_CHECK_EQUAL(tst_obj.ts_flux_diffuse_colwf().extent(2), lid_pars.max_szangles);
  BOOST_CHECK_EQUAL(tst_obj.ts_flux_diffuse_colwf().extent(3), lid_pars.maxstokes);
  BOOST_CHECK_EQUAL(tst_obj.ts_flux_diffuse_colwf().extent(4), lid_pars.max_directions);
  BOOST_CHECK_EQUAL(tst_obj.ts_dnmeanst_direct_colwf().extent(0), lid_pars.max_atmoswfs);
  BOOST_CHECK_EQUAL(tst_obj.ts_dnmeanst_direct_colwf().extent(1), lid_pars.max_user_levels);
  BOOST_CHECK_EQUAL(tst_obj.ts_dnmeanst_direct_colwf().extent(2), lid_pars.max_szangles);
  BOOST_CHECK_EQUAL(tst_obj.ts_dnmeanst_direct_colwf().extent(3), lid_pars.maxstokes);
  BOOST_CHECK_EQUAL(tst_obj.ts_dnflux_direct_colwf().extent(0), lid_pars.max_atmoswfs);
  BOOST_CHECK_EQUAL(tst_obj.ts_dnflux_direct_colwf().extent(1), lid_pars.max_user_levels);
  BOOST_CHECK_EQUAL(tst_obj.ts_dnflux_direct_colwf().extent(2), lid_pars.max_szangles);
  BOOST_CHECK_EQUAL(tst_obj.ts_dnflux_direct_colwf().extent(3), lid_pars.maxstokes);
  BOOST_CHECK_EQUAL(tst_obj.ts_albmed_user_colwf().extent(0), lid_pars.maxstokes);
  BOOST_CHECK_EQUAL(tst_obj.ts_albmed_user_colwf().extent(1), lid_pars.max_user_vzangles);
  BOOST_CHECK_EQUAL(tst_obj.ts_albmed_user_colwf().extent(2), lid_pars.max_atmoswfs);
  BOOST_CHECK_EQUAL(tst_obj.ts_trnmed_user_colwf().extent(0), lid_pars.maxstokes);
  BOOST_CHECK_EQUAL(tst_obj.ts_trnmed_user_colwf().extent(1), lid_pars.max_user_vzangles);
  BOOST_CHECK_EQUAL(tst_obj.ts_trnmed_user_colwf().extent(2), lid_pars.max_atmoswfs);
  BOOST_CHECK_EQUAL(tst_obj.ts_albmed_fluxes_colwf().extent(0), lid_pars.maxstokes);
  BOOST_CHECK_EQUAL(tst_obj.ts_albmed_fluxes_colwf().extent(1), lid_pars.2);
  BOOST_CHECK_EQUAL(tst_obj.ts_albmed_fluxes_colwf().extent(2), lid_pars.max_atmoswfs);
  BOOST_CHECK_EQUAL(tst_obj.ts_trnmed_fluxes_colwf().extent(0), lid_pars.maxstokes);
  BOOST_CHECK_EQUAL(tst_obj.ts_trnmed_fluxes_colwf().extent(1), lid_pars.2);
  BOOST_CHECK_EQUAL(tst_obj.ts_trnmed_fluxes_colwf().extent(2), lid_pars.max_atmoswfs);
  BOOST_CHECK_EQUAL(tst_obj.ts_transbeam_colwf().extent(0), lid_pars.maxstokes);
  BOOST_CHECK_EQUAL(tst_obj.ts_transbeam_colwf().extent(1), lid_pars.max_szangles);
  BOOST_CHECK_EQUAL(tst_obj.ts_transbeam_colwf().extent(2), lid_pars.max_atmoswfs);
  BOOST_CHECK_EQUAL(tst_obj.ts_planetary_transterm_colwf().extent(0), lid_pars.maxstokes);
  BOOST_CHECK_EQUAL(tst_obj.ts_planetary_transterm_colwf().extent(1), lid_pars.max_geometries);
  BOOST_CHECK_EQUAL(tst_obj.ts_planetary_transterm_colwf().extent(2), lid_pars.max_atmoswfs);
  BOOST_CHECK_EQUAL(tst_obj.ts_planetary_sbterm_colwf().extent(0), lid_pars.max_atmoswfs);
  BOOST_CHECK_EQUAL(tst_obj.ts_lc_lostrans().extent(0), lid_pars.max_atmoswfs);
  BOOST_CHECK_EQUAL(tst_obj.ts_lc_lostrans().extent(1), lid_pars.max_szangles);
  BOOST_CHECK_EQUAL(tst_obj.ts_lc_lostrans().extent(2), lid_pars.maxlayers);
  BOOST_CHECK_EQUAL(tst_obj.ts_lc_layer_mssts().extent(0), lid_pars.max_atmoswfs);
  BOOST_CHECK_EQUAL(tst_obj.ts_lc_layer_mssts().extent(1), lid_pars.max_szangles);
  BOOST_CHECK_EQUAL(tst_obj.ts_lc_layer_mssts().extent(2), lid_pars.maxstokes);
  BOOST_CHECK_EQUAL(tst_obj.ts_lc_layer_mssts().extent(3), lid_pars.maxlayers);
  BOOST_CHECK_EQUAL(tst_obj.ts_lc_surf_mssts().extent(0), lid_pars.max_atmoswfs);
  BOOST_CHECK_EQUAL(tst_obj.ts_lc_surf_mssts().extent(1), lid_pars.max_szangles);
  BOOST_CHECK_EQUAL(tst_obj.ts_lc_surf_mssts().extent(2), lid_pars.maxstokes);
  

  // Test initialization
  blitz::Array<double, 5> ts_columnwf_exp(tst_obj.ts_columnwf().shape());
  ts_columnwf_exp = 0;
  BOOST_CHECK_MATRIX_CLOSE_TOL(tst_obj.ts_columnwf(), ts_columnwf_exp, 1e-10);
  blitz::Array<double, 5> ts_meanst_diffuse_colwf_exp(tst_obj.ts_meanst_diffuse_colwf().shape());
  ts_meanst_diffuse_colwf_exp = 0;
  BOOST_CHECK_MATRIX_CLOSE_TOL(tst_obj.ts_meanst_diffuse_colwf(), ts_meanst_diffuse_colwf_exp, 1e-10);
  blitz::Array<double, 5> ts_flux_diffuse_colwf_exp(tst_obj.ts_flux_diffuse_colwf().shape());
  ts_flux_diffuse_colwf_exp = 0;
  BOOST_CHECK_MATRIX_CLOSE_TOL(tst_obj.ts_flux_diffuse_colwf(), ts_flux_diffuse_colwf_exp, 1e-10);
  blitz::Array<double, 4> ts_dnmeanst_direct_colwf_exp(tst_obj.ts_dnmeanst_direct_colwf().shape());
  ts_dnmeanst_direct_colwf_exp = 0;
  BOOST_CHECK_MATRIX_CLOSE_TOL(tst_obj.ts_dnmeanst_direct_colwf(), ts_dnmeanst_direct_colwf_exp, 1e-10);
  blitz::Array<double, 4> ts_dnflux_direct_colwf_exp(tst_obj.ts_dnflux_direct_colwf().shape());
  ts_dnflux_direct_colwf_exp = 0;
  BOOST_CHECK_MATRIX_CLOSE_TOL(tst_obj.ts_dnflux_direct_colwf(), ts_dnflux_direct_colwf_exp, 1e-10);
  blitz::Array<double, 3> ts_albmed_user_colwf_exp(tst_obj.ts_albmed_user_colwf().shape());
  ts_albmed_user_colwf_exp = 0;
  BOOST_CHECK_MATRIX_CLOSE_TOL(tst_obj.ts_albmed_user_colwf(), ts_albmed_user_colwf_exp, 1e-10);
  blitz::Array<double, 3> ts_trnmed_user_colwf_exp(tst_obj.ts_trnmed_user_colwf().shape());
  ts_trnmed_user_colwf_exp = 0;
  BOOST_CHECK_MATRIX_CLOSE_TOL(tst_obj.ts_trnmed_user_colwf(), ts_trnmed_user_colwf_exp, 1e-10);
  blitz::Array<double, 3> ts_albmed_fluxes_colwf_exp(tst_obj.ts_albmed_fluxes_colwf().shape());
  ts_albmed_fluxes_colwf_exp = 0;
  BOOST_CHECK_MATRIX_CLOSE_TOL(tst_obj.ts_albmed_fluxes_colwf(), ts_albmed_fluxes_colwf_exp, 1e-10);
  blitz::Array<double, 3> ts_trnmed_fluxes_colwf_exp(tst_obj.ts_trnmed_fluxes_colwf().shape());
  ts_trnmed_fluxes_colwf_exp = 0;
  BOOST_CHECK_MATRIX_CLOSE_TOL(tst_obj.ts_trnmed_fluxes_colwf(), ts_trnmed_fluxes_colwf_exp, 1e-10);
  blitz::Array<double, 3> ts_transbeam_colwf_exp(tst_obj.ts_transbeam_colwf().shape());
  ts_transbeam_colwf_exp = 0;
  BOOST_CHECK_MATRIX_CLOSE_TOL(tst_obj.ts_transbeam_colwf(), ts_transbeam_colwf_exp, 1e-10);
  blitz::Array<double, 3> ts_planetary_transterm_colwf_exp(tst_obj.ts_planetary_transterm_colwf().shape());
  ts_planetary_transterm_colwf_exp = 0;
  BOOST_CHECK_MATRIX_CLOSE_TOL(tst_obj.ts_planetary_transterm_colwf(), ts_planetary_transterm_colwf_exp, 1e-10);
  blitz::Array<double, 1> ts_planetary_sbterm_colwf_exp(tst_obj.ts_planetary_sbterm_colwf().shape());
  ts_planetary_sbterm_colwf_exp = 0;
  BOOST_CHECK_MATRIX_CLOSE_TOL(tst_obj.ts_planetary_sbterm_colwf(), ts_planetary_sbterm_colwf_exp, 1e-10);
  blitz::Array<double, 3> ts_lc_lostrans_exp(tst_obj.ts_lc_lostrans().shape());
  ts_lc_lostrans_exp = 0;
  BOOST_CHECK_MATRIX_CLOSE_TOL(tst_obj.ts_lc_lostrans(), ts_lc_lostrans_exp, 1e-10);
  blitz::Array<double, 4> ts_lc_layer_mssts_exp(tst_obj.ts_lc_layer_mssts().shape());
  ts_lc_layer_mssts_exp = 0;
  BOOST_CHECK_MATRIX_CLOSE_TOL(tst_obj.ts_lc_layer_mssts(), ts_lc_layer_mssts_exp, 1e-10);
  blitz::Array<double, 3> ts_lc_surf_mssts_exp(tst_obj.ts_lc_surf_mssts().shape());
  ts_lc_surf_mssts_exp = 0;
  BOOST_CHECK_MATRIX_CLOSE_TOL(tst_obj.ts_lc_surf_mssts(), ts_lc_surf_mssts_exp, 1e-10);
  
}

BOOST_AUTO_TEST_CASE(vlidort_linprof)
{
  // Used for checking dimensions
  VLidort_Pars lid_pars = VLidort_Pars::instance();

  // Test constructor
  VLidort_Linprof tst_obj = VLidort_Linprof();

  // Test memory space sizes
  tst_obj.check_byte_sizes(); // throws exception on mismatch

  // Test variable shapes
  BOOST_CHECK_EQUAL(tst_obj.ts_profilewf().extent(0), lid_pars.max_atmoswfs);
  BOOST_CHECK_EQUAL(tst_obj.ts_profilewf().extent(1), lid_pars.maxlayers);
  BOOST_CHECK_EQUAL(tst_obj.ts_profilewf().extent(2), lid_pars.max_user_levels);
  BOOST_CHECK_EQUAL(tst_obj.ts_profilewf().extent(3), lid_pars.max_geometries);
  BOOST_CHECK_EQUAL(tst_obj.ts_profilewf().extent(4), lid_pars.maxstokes);
  BOOST_CHECK_EQUAL(tst_obj.ts_profilewf().extent(5), lid_pars.max_directions);
  BOOST_CHECK_EQUAL(tst_obj.ts_meanst_diffuse_profwf().extent(0), lid_pars.max_atmoswfs);
  BOOST_CHECK_EQUAL(tst_obj.ts_meanst_diffuse_profwf().extent(1), lid_pars.maxlayers);
  BOOST_CHECK_EQUAL(tst_obj.ts_meanst_diffuse_profwf().extent(2), lid_pars.max_user_levels);
  BOOST_CHECK_EQUAL(tst_obj.ts_meanst_diffuse_profwf().extent(3), lid_pars.max_szangles);
  BOOST_CHECK_EQUAL(tst_obj.ts_meanst_diffuse_profwf().extent(4), lid_pars.maxstokes);
  BOOST_CHECK_EQUAL(tst_obj.ts_meanst_diffuse_profwf().extent(5), lid_pars.max_directions);
  BOOST_CHECK_EQUAL(tst_obj.ts_flux_diffuse_profwf().extent(0), lid_pars.max_atmoswfs);
  BOOST_CHECK_EQUAL(tst_obj.ts_flux_diffuse_profwf().extent(1), lid_pars.maxlayers);
  BOOST_CHECK_EQUAL(tst_obj.ts_flux_diffuse_profwf().extent(2), lid_pars.max_user_levels);
  BOOST_CHECK_EQUAL(tst_obj.ts_flux_diffuse_profwf().extent(3), lid_pars.max_szangles);
  BOOST_CHECK_EQUAL(tst_obj.ts_flux_diffuse_profwf().extent(4), lid_pars.maxstokes);
  BOOST_CHECK_EQUAL(tst_obj.ts_flux_diffuse_profwf().extent(5), lid_pars.max_directions);
  BOOST_CHECK_EQUAL(tst_obj.ts_dnmeanst_direct_profwf().extent(0), lid_pars.max_atmoswfs);
  BOOST_CHECK_EQUAL(tst_obj.ts_dnmeanst_direct_profwf().extent(1), lid_pars.maxlayers);
  BOOST_CHECK_EQUAL(tst_obj.ts_dnmeanst_direct_profwf().extent(2), lid_pars.max_user_levels);
  BOOST_CHECK_EQUAL(tst_obj.ts_dnmeanst_direct_profwf().extent(3), lid_pars.max_szangles);
  BOOST_CHECK_EQUAL(tst_obj.ts_dnmeanst_direct_profwf().extent(4), lid_pars.maxstokes);
  BOOST_CHECK_EQUAL(tst_obj.ts_dnflux_direct_profwf().extent(0), lid_pars.max_atmoswfs);
  BOOST_CHECK_EQUAL(tst_obj.ts_dnflux_direct_profwf().extent(1), lid_pars.maxlayers);
  BOOST_CHECK_EQUAL(tst_obj.ts_dnflux_direct_profwf().extent(2), lid_pars.max_user_levels);
  BOOST_CHECK_EQUAL(tst_obj.ts_dnflux_direct_profwf().extent(3), lid_pars.max_szangles);
  BOOST_CHECK_EQUAL(tst_obj.ts_dnflux_direct_profwf().extent(4), lid_pars.maxstokes);
  BOOST_CHECK_EQUAL(tst_obj.ts_albmed_user_profwf().extent(0), lid_pars.maxstokes);
  BOOST_CHECK_EQUAL(tst_obj.ts_albmed_user_profwf().extent(1), lid_pars.max_user_vzangles);
  BOOST_CHECK_EQUAL(tst_obj.ts_albmed_user_profwf().extent(2), lid_pars.maxlayers);
  BOOST_CHECK_EQUAL(tst_obj.ts_albmed_user_profwf().extent(3), lid_pars.max_atmoswfs);
  BOOST_CHECK_EQUAL(tst_obj.ts_trnmed_user_profwf().extent(0), lid_pars.maxstokes);
  BOOST_CHECK_EQUAL(tst_obj.ts_trnmed_user_profwf().extent(1), lid_pars.max_user_vzangles);
  BOOST_CHECK_EQUAL(tst_obj.ts_trnmed_user_profwf().extent(2), lid_pars.maxlayers);
  BOOST_CHECK_EQUAL(tst_obj.ts_trnmed_user_profwf().extent(3), lid_pars.max_atmoswfs);
  BOOST_CHECK_EQUAL(tst_obj.ts_albmed_fluxes_profwf().extent(0), lid_pars.maxstokes);
  BOOST_CHECK_EQUAL(tst_obj.ts_albmed_fluxes_profwf().extent(1), lid_pars.2);
  BOOST_CHECK_EQUAL(tst_obj.ts_albmed_fluxes_profwf().extent(2), lid_pars.maxlayers);
  BOOST_CHECK_EQUAL(tst_obj.ts_albmed_fluxes_profwf().extent(3), lid_pars.max_atmoswfs);
  BOOST_CHECK_EQUAL(tst_obj.ts_trnmed_fluxes_profwf().extent(0), lid_pars.maxstokes);
  BOOST_CHECK_EQUAL(tst_obj.ts_trnmed_fluxes_profwf().extent(1), lid_pars.2);
  BOOST_CHECK_EQUAL(tst_obj.ts_trnmed_fluxes_profwf().extent(2), lid_pars.maxlayers);
  BOOST_CHECK_EQUAL(tst_obj.ts_trnmed_fluxes_profwf().extent(3), lid_pars.max_atmoswfs);
  BOOST_CHECK_EQUAL(tst_obj.ts_transbeam_profwf().extent(0), lid_pars.maxstokes);
  BOOST_CHECK_EQUAL(tst_obj.ts_transbeam_profwf().extent(1), lid_pars.max_szangles);
  BOOST_CHECK_EQUAL(tst_obj.ts_transbeam_profwf().extent(2), lid_pars.maxlayers);
  BOOST_CHECK_EQUAL(tst_obj.ts_transbeam_profwf().extent(3), lid_pars.max_atmoswfs);
  BOOST_CHECK_EQUAL(tst_obj.ts_planetary_transterm_profwf().extent(0), lid_pars.maxstokes);
  BOOST_CHECK_EQUAL(tst_obj.ts_planetary_transterm_profwf().extent(1), lid_pars.max_geometries);
  BOOST_CHECK_EQUAL(tst_obj.ts_planetary_transterm_profwf().extent(2), lid_pars.maxlayers);
  BOOST_CHECK_EQUAL(tst_obj.ts_planetary_transterm_profwf().extent(3), lid_pars.max_atmoswfs);
  BOOST_CHECK_EQUAL(tst_obj.ts_planetary_sbterm_profwf().extent(0), lid_pars.maxlayers);
  BOOST_CHECK_EQUAL(tst_obj.ts_planetary_sbterm_profwf().extent(1), lid_pars.max_atmoswfs);
  BOOST_CHECK_EQUAL(tst_obj.ts_lp_lostrans().extent(0), lid_pars.max_atmoswfs);
  BOOST_CHECK_EQUAL(tst_obj.ts_lp_lostrans().extent(1), lid_pars.max_szangles);
  BOOST_CHECK_EQUAL(tst_obj.ts_lp_lostrans().extent(2), lid_pars.maxlayers);
  BOOST_CHECK_EQUAL(tst_obj.ts_lp_layer_mssts().extent(0), lid_pars.max_atmoswfs);
  BOOST_CHECK_EQUAL(tst_obj.ts_lp_layer_mssts().extent(1), lid_pars.maxlayers);
  BOOST_CHECK_EQUAL(tst_obj.ts_lp_layer_mssts().extent(2), lid_pars.max_szangles);
  BOOST_CHECK_EQUAL(tst_obj.ts_lp_layer_mssts().extent(3), lid_pars.maxstokes);
  BOOST_CHECK_EQUAL(tst_obj.ts_lp_layer_mssts().extent(4), lid_pars.maxlayers);
  BOOST_CHECK_EQUAL(tst_obj.ts_lp_surf_mssts().extent(0), lid_pars.max_atmoswfs);
  BOOST_CHECK_EQUAL(tst_obj.ts_lp_surf_mssts().extent(1), lid_pars.maxlayers);
  BOOST_CHECK_EQUAL(tst_obj.ts_lp_surf_mssts().extent(2), lid_pars.max_szangles);
  BOOST_CHECK_EQUAL(tst_obj.ts_lp_surf_mssts().extent(3), lid_pars.maxstokes);
  

  // Test initialization
  blitz::Array<double, 6> ts_profilewf_exp(tst_obj.ts_profilewf().shape());
  ts_profilewf_exp = 0;
  BOOST_CHECK_MATRIX_CLOSE_TOL(tst_obj.ts_profilewf(), ts_profilewf_exp, 1e-10);
  blitz::Array<double, 6> ts_meanst_diffuse_profwf_exp(tst_obj.ts_meanst_diffuse_profwf().shape());
  ts_meanst_diffuse_profwf_exp = 0;
  BOOST_CHECK_MATRIX_CLOSE_TOL(tst_obj.ts_meanst_diffuse_profwf(), ts_meanst_diffuse_profwf_exp, 1e-10);
  blitz::Array<double, 6> ts_flux_diffuse_profwf_exp(tst_obj.ts_flux_diffuse_profwf().shape());
  ts_flux_diffuse_profwf_exp = 0;
  BOOST_CHECK_MATRIX_CLOSE_TOL(tst_obj.ts_flux_diffuse_profwf(), ts_flux_diffuse_profwf_exp, 1e-10);
  blitz::Array<double, 5> ts_dnmeanst_direct_profwf_exp(tst_obj.ts_dnmeanst_direct_profwf().shape());
  ts_dnmeanst_direct_profwf_exp = 0;
  BOOST_CHECK_MATRIX_CLOSE_TOL(tst_obj.ts_dnmeanst_direct_profwf(), ts_dnmeanst_direct_profwf_exp, 1e-10);
  blitz::Array<double, 5> ts_dnflux_direct_profwf_exp(tst_obj.ts_dnflux_direct_profwf().shape());
  ts_dnflux_direct_profwf_exp = 0;
  BOOST_CHECK_MATRIX_CLOSE_TOL(tst_obj.ts_dnflux_direct_profwf(), ts_dnflux_direct_profwf_exp, 1e-10);
  blitz::Array<double, 4> ts_albmed_user_profwf_exp(tst_obj.ts_albmed_user_profwf().shape());
  ts_albmed_user_profwf_exp = 0;
  BOOST_CHECK_MATRIX_CLOSE_TOL(tst_obj.ts_albmed_user_profwf(), ts_albmed_user_profwf_exp, 1e-10);
  blitz::Array<double, 4> ts_trnmed_user_profwf_exp(tst_obj.ts_trnmed_user_profwf().shape());
  ts_trnmed_user_profwf_exp = 0;
  BOOST_CHECK_MATRIX_CLOSE_TOL(tst_obj.ts_trnmed_user_profwf(), ts_trnmed_user_profwf_exp, 1e-10);
  blitz::Array<double, 4> ts_albmed_fluxes_profwf_exp(tst_obj.ts_albmed_fluxes_profwf().shape());
  ts_albmed_fluxes_profwf_exp = 0;
  BOOST_CHECK_MATRIX_CLOSE_TOL(tst_obj.ts_albmed_fluxes_profwf(), ts_albmed_fluxes_profwf_exp, 1e-10);
  blitz::Array<double, 4> ts_trnmed_fluxes_profwf_exp(tst_obj.ts_trnmed_fluxes_profwf().shape());
  ts_trnmed_fluxes_profwf_exp = 0;
  BOOST_CHECK_MATRIX_CLOSE_TOL(tst_obj.ts_trnmed_fluxes_profwf(), ts_trnmed_fluxes_profwf_exp, 1e-10);
  blitz::Array<double, 4> ts_transbeam_profwf_exp(tst_obj.ts_transbeam_profwf().shape());
  ts_transbeam_profwf_exp = 0;
  BOOST_CHECK_MATRIX_CLOSE_TOL(tst_obj.ts_transbeam_profwf(), ts_transbeam_profwf_exp, 1e-10);
  blitz::Array<double, 4> ts_planetary_transterm_profwf_exp(tst_obj.ts_planetary_transterm_profwf().shape());
  ts_planetary_transterm_profwf_exp = 0;
  BOOST_CHECK_MATRIX_CLOSE_TOL(tst_obj.ts_planetary_transterm_profwf(), ts_planetary_transterm_profwf_exp, 1e-10);
  blitz::Array<double, 2> ts_planetary_sbterm_profwf_exp(tst_obj.ts_planetary_sbterm_profwf().shape());
  ts_planetary_sbterm_profwf_exp = 0;
  BOOST_CHECK_MATRIX_CLOSE_TOL(tst_obj.ts_planetary_sbterm_profwf(), ts_planetary_sbterm_profwf_exp, 1e-10);
  blitz::Array<double, 3> ts_lp_lostrans_exp(tst_obj.ts_lp_lostrans().shape());
  ts_lp_lostrans_exp = 0;
  BOOST_CHECK_MATRIX_CLOSE_TOL(tst_obj.ts_lp_lostrans(), ts_lp_lostrans_exp, 1e-10);
  blitz::Array<double, 5> ts_lp_layer_mssts_exp(tst_obj.ts_lp_layer_mssts().shape());
  ts_lp_layer_mssts_exp = 0;
  BOOST_CHECK_MATRIX_CLOSE_TOL(tst_obj.ts_lp_layer_mssts(), ts_lp_layer_mssts_exp, 1e-10);
  blitz::Array<double, 4> ts_lp_surf_mssts_exp(tst_obj.ts_lp_surf_mssts().shape());
  ts_lp_surf_mssts_exp = 0;
  BOOST_CHECK_MATRIX_CLOSE_TOL(tst_obj.ts_lp_surf_mssts(), ts_lp_surf_mssts_exp, 1e-10);
  
}

BOOST_AUTO_TEST_CASE(vlidort_linatmos)
{
  // Used for checking dimensions
  VLidort_Pars lid_pars = VLidort_Pars::instance();

  // Test constructor
  VLidort_Linatmos tst_obj = VLidort_Linatmos();

  // Test memory space sizes
  tst_obj.check_byte_sizes(); // throws exception on mismatch

  // Test variable shapes
  BOOST_CHECK_EQUAL(tst_obj.ts_abbwfs_jacobians().extent(0), lid_pars.max_user_levels);
  BOOST_CHECK_EQUAL(tst_obj.ts_abbwfs_jacobians().extent(1), lid_pars.max_user_vzangles);
  BOOST_CHECK_EQUAL(tst_obj.ts_abbwfs_jacobians().extent(2), lid_pars.maxlayers+1);
  BOOST_CHECK_EQUAL(tst_obj.ts_abbwfs_jacobians().extent(3), lid_pars.maxstokes);
  BOOST_CHECK_EQUAL(tst_obj.ts_abbwfs_jacobians().extent(4), lid_pars.max_directions);
  BOOST_CHECK_EQUAL(tst_obj.ts_abbwfs_fluxes().extent(0), lid_pars.max_user_levels);
  BOOST_CHECK_EQUAL(tst_obj.ts_abbwfs_fluxes().extent(1), lid_pars.2);
  BOOST_CHECK_EQUAL(tst_obj.ts_abbwfs_fluxes().extent(2), lid_pars.maxlayers+1);
  BOOST_CHECK_EQUAL(tst_obj.ts_abbwfs_fluxes().extent(3), lid_pars.maxstokes);
  BOOST_CHECK_EQUAL(tst_obj.ts_abbwfs_fluxes().extent(4), lid_pars.max_directions);
  

  // Test initialization
  blitz::Array<double, 5> ts_abbwfs_jacobians_exp(tst_obj.ts_abbwfs_jacobians().shape());
  ts_abbwfs_jacobians_exp = 0;
  BOOST_CHECK_MATRIX_CLOSE_TOL(tst_obj.ts_abbwfs_jacobians(), ts_abbwfs_jacobians_exp, 1e-10);
  blitz::Array<double, 5> ts_abbwfs_fluxes_exp(tst_obj.ts_abbwfs_fluxes().shape());
  ts_abbwfs_fluxes_exp = 0;
  BOOST_CHECK_MATRIX_CLOSE_TOL(tst_obj.ts_abbwfs_fluxes(), ts_abbwfs_fluxes_exp, 1e-10);
  
}

BOOST_AUTO_TEST_CASE(vlidort_linsurf)
{
  // Used for checking dimensions
  VLidort_Pars lid_pars = VLidort_Pars::instance();

  // Test constructor
  VLidort_Linsurf tst_obj = VLidort_Linsurf();

  // Test memory space sizes
  tst_obj.check_byte_sizes(); // throws exception on mismatch

  // Test variable shapes
  BOOST_CHECK_EQUAL(tst_obj.ts_surfacewf().extent(0), lid_pars.max_surfacewfs);
  BOOST_CHECK_EQUAL(tst_obj.ts_surfacewf().extent(1), lid_pars.max_user_levels);
  BOOST_CHECK_EQUAL(tst_obj.ts_surfacewf().extent(2), lid_pars.max_geometries);
  BOOST_CHECK_EQUAL(tst_obj.ts_surfacewf().extent(3), lid_pars.maxstokes);
  BOOST_CHECK_EQUAL(tst_obj.ts_surfacewf().extent(4), lid_pars.max_directions);
  BOOST_CHECK_EQUAL(tst_obj.ts_meanst_diffuse_surfwf().extent(0), lid_pars.max_surfacewfs);
  BOOST_CHECK_EQUAL(tst_obj.ts_meanst_diffuse_surfwf().extent(1), lid_pars.max_user_levels);
  BOOST_CHECK_EQUAL(tst_obj.ts_meanst_diffuse_surfwf().extent(2), lid_pars.max_szangles);
  BOOST_CHECK_EQUAL(tst_obj.ts_meanst_diffuse_surfwf().extent(3), lid_pars.maxstokes);
  BOOST_CHECK_EQUAL(tst_obj.ts_meanst_diffuse_surfwf().extent(4), lid_pars.max_directions);
  BOOST_CHECK_EQUAL(tst_obj.ts_flux_diffuse_surfwf().extent(0), lid_pars.max_surfacewfs);
  BOOST_CHECK_EQUAL(tst_obj.ts_flux_diffuse_surfwf().extent(1), lid_pars.max_user_levels);
  BOOST_CHECK_EQUAL(tst_obj.ts_flux_diffuse_surfwf().extent(2), lid_pars.max_szangles);
  BOOST_CHECK_EQUAL(tst_obj.ts_flux_diffuse_surfwf().extent(3), lid_pars.maxstokes);
  BOOST_CHECK_EQUAL(tst_obj.ts_flux_diffuse_surfwf().extent(4), lid_pars.max_directions);
  BOOST_CHECK_EQUAL(tst_obj.ts_ls_layer_mssts().extent(0), lid_pars.max_surfacewfs);
  BOOST_CHECK_EQUAL(tst_obj.ts_ls_layer_mssts().extent(1), lid_pars.max_szangles);
  BOOST_CHECK_EQUAL(tst_obj.ts_ls_layer_mssts().extent(2), lid_pars.maxstokes);
  BOOST_CHECK_EQUAL(tst_obj.ts_ls_layer_mssts().extent(3), lid_pars.maxlayers);
  BOOST_CHECK_EQUAL(tst_obj.ts_ls_surf_mssts().extent(0), lid_pars.max_surfacewfs);
  BOOST_CHECK_EQUAL(tst_obj.ts_ls_surf_mssts().extent(1), lid_pars.max_szangles);
  BOOST_CHECK_EQUAL(tst_obj.ts_ls_surf_mssts().extent(2), lid_pars.maxstokes);
  BOOST_CHECK_EQUAL(tst_obj.ts_sbbwfs_jacobians().extent(0), lid_pars.max_user_levels);
  BOOST_CHECK_EQUAL(tst_obj.ts_sbbwfs_jacobians().extent(1), lid_pars.max_user_vzangles);
  BOOST_CHECK_EQUAL(tst_obj.ts_sbbwfs_jacobians().extent(2), lid_pars.maxstokes);
  BOOST_CHECK_EQUAL(tst_obj.ts_sbbwfs_jacobians().extent(3), lid_pars.max_directions);
  BOOST_CHECK_EQUAL(tst_obj.ts_sbbwfs_fluxes().extent(0), lid_pars.max_user_levels);
  BOOST_CHECK_EQUAL(tst_obj.ts_sbbwfs_fluxes().extent(1), lid_pars.2);
  BOOST_CHECK_EQUAL(tst_obj.ts_sbbwfs_fluxes().extent(2), lid_pars.maxstokes);
  BOOST_CHECK_EQUAL(tst_obj.ts_sbbwfs_fluxes().extent(3), lid_pars.max_directions);
  

  // Test initialization
  blitz::Array<double, 5> ts_surfacewf_exp(tst_obj.ts_surfacewf().shape());
  ts_surfacewf_exp = 0;
  BOOST_CHECK_MATRIX_CLOSE_TOL(tst_obj.ts_surfacewf(), ts_surfacewf_exp, 1e-10);
  blitz::Array<double, 5> ts_meanst_diffuse_surfwf_exp(tst_obj.ts_meanst_diffuse_surfwf().shape());
  ts_meanst_diffuse_surfwf_exp = 0;
  BOOST_CHECK_MATRIX_CLOSE_TOL(tst_obj.ts_meanst_diffuse_surfwf(), ts_meanst_diffuse_surfwf_exp, 1e-10);
  blitz::Array<double, 5> ts_flux_diffuse_surfwf_exp(tst_obj.ts_flux_diffuse_surfwf().shape());
  ts_flux_diffuse_surfwf_exp = 0;
  BOOST_CHECK_MATRIX_CLOSE_TOL(tst_obj.ts_flux_diffuse_surfwf(), ts_flux_diffuse_surfwf_exp, 1e-10);
  blitz::Array<double, 4> ts_ls_layer_mssts_exp(tst_obj.ts_ls_layer_mssts().shape());
  ts_ls_layer_mssts_exp = 0;
  BOOST_CHECK_MATRIX_CLOSE_TOL(tst_obj.ts_ls_layer_mssts(), ts_ls_layer_mssts_exp, 1e-10);
  blitz::Array<double, 3> ts_ls_surf_mssts_exp(tst_obj.ts_ls_surf_mssts().shape());
  ts_ls_surf_mssts_exp = 0;
  BOOST_CHECK_MATRIX_CLOSE_TOL(tst_obj.ts_ls_surf_mssts(), ts_ls_surf_mssts_exp, 1e-10);
  blitz::Array<double, 4> ts_sbbwfs_jacobians_exp(tst_obj.ts_sbbwfs_jacobians().shape());
  ts_sbbwfs_jacobians_exp = 0;
  BOOST_CHECK_MATRIX_CLOSE_TOL(tst_obj.ts_sbbwfs_jacobians(), ts_sbbwfs_jacobians_exp, 1e-10);
  blitz::Array<double, 4> ts_sbbwfs_fluxes_exp(tst_obj.ts_sbbwfs_fluxes().shape());
  ts_sbbwfs_fluxes_exp = 0;
  BOOST_CHECK_MATRIX_CLOSE_TOL(tst_obj.ts_sbbwfs_fluxes(), ts_sbbwfs_fluxes_exp, 1e-10);
  
}

BOOST_AUTO_TEST_CASE(vlidort_linoutputs)
{
  // Used for checking dimensions
  VLidort_Pars lid_pars = VLidort_Pars::instance();

  // Test constructor
  VLidort_Linoutputs tst_obj = VLidort_Linoutputs();

  // Test memory space sizes
  tst_obj.check_byte_sizes(); // throws exception on mismatch

  // Test variable shapes
  

  // Test initialization
  
}

BOOST_AUTO_TEST_CASE(vlidort_linsup_brdf)
{
  // Used for checking dimensions
  VLidort_Pars lid_pars = VLidort_Pars::instance();

  // Test constructor
  VLidort_Linsup_Brdf tst_obj = VLidort_Linsup_Brdf();

  // Test memory space sizes
  tst_obj.check_byte_sizes(); // throws exception on mismatch

  // Test variable shapes
  BOOST_CHECK_EQUAL(tst_obj.ts_ls_exactdb_brdfunc().extent(0), lid_pars.max_surfacewfs);
  BOOST_CHECK_EQUAL(tst_obj.ts_ls_exactdb_brdfunc().extent(1), lid_pars.maxstokes_sq);
  BOOST_CHECK_EQUAL(tst_obj.ts_ls_exactdb_brdfunc().extent(2), lid_pars.max_user_streams);
  BOOST_CHECK_EQUAL(tst_obj.ts_ls_exactdb_brdfunc().extent(3), lid_pars.max_user_relazms);
  BOOST_CHECK_EQUAL(tst_obj.ts_ls_exactdb_brdfunc().extent(4), lid_pars.maxbeams);
  BOOST_CHECK_EQUAL(tst_obj.ts_ls_brdf_f_0().extent(0), lid_pars.max_surfacewfs);
  BOOST_CHECK_EQUAL(tst_obj.ts_ls_brdf_f_0().extent(1), lid_pars.maxmoments+1);
  BOOST_CHECK_EQUAL(tst_obj.ts_ls_brdf_f_0().extent(2), lid_pars.maxstokes_sq);
  BOOST_CHECK_EQUAL(tst_obj.ts_ls_brdf_f_0().extent(3), lid_pars.maxstreams);
  BOOST_CHECK_EQUAL(tst_obj.ts_ls_brdf_f_0().extent(4), lid_pars.maxbeams);
  BOOST_CHECK_EQUAL(tst_obj.ts_ls_brdf_f().extent(0), lid_pars.max_surfacewfs);
  BOOST_CHECK_EQUAL(tst_obj.ts_ls_brdf_f().extent(1), lid_pars.maxmoments+1);
  BOOST_CHECK_EQUAL(tst_obj.ts_ls_brdf_f().extent(2), lid_pars.maxstokes_sq);
  BOOST_CHECK_EQUAL(tst_obj.ts_ls_brdf_f().extent(3), lid_pars.maxstreams);
  BOOST_CHECK_EQUAL(tst_obj.ts_ls_brdf_f().extent(4), lid_pars.maxstreams);
  BOOST_CHECK_EQUAL(tst_obj.ts_ls_user_brdf_f_0().extent(0), lid_pars.max_surfacewfs);
  BOOST_CHECK_EQUAL(tst_obj.ts_ls_user_brdf_f_0().extent(1), lid_pars.maxmoments+1);
  BOOST_CHECK_EQUAL(tst_obj.ts_ls_user_brdf_f_0().extent(2), lid_pars.maxstokes_sq);
  BOOST_CHECK_EQUAL(tst_obj.ts_ls_user_brdf_f_0().extent(3), lid_pars.max_user_streams);
  BOOST_CHECK_EQUAL(tst_obj.ts_ls_user_brdf_f_0().extent(4), lid_pars.maxbeams);
  BOOST_CHECK_EQUAL(tst_obj.ts_ls_user_brdf_f().extent(0), lid_pars.max_surfacewfs);
  BOOST_CHECK_EQUAL(tst_obj.ts_ls_user_brdf_f().extent(1), lid_pars.maxmoments+1);
  BOOST_CHECK_EQUAL(tst_obj.ts_ls_user_brdf_f().extent(2), lid_pars.maxstokes_sq);
  BOOST_CHECK_EQUAL(tst_obj.ts_ls_user_brdf_f().extent(3), lid_pars.max_user_streams);
  BOOST_CHECK_EQUAL(tst_obj.ts_ls_user_brdf_f().extent(4), lid_pars.maxstreams);
  BOOST_CHECK_EQUAL(tst_obj.ts_ls_emissivity().extent(0), lid_pars.max_surfacewfs);
  BOOST_CHECK_EQUAL(tst_obj.ts_ls_emissivity().extent(1), lid_pars.maxstokes);
  BOOST_CHECK_EQUAL(tst_obj.ts_ls_emissivity().extent(2), lid_pars.maxstreams);
  BOOST_CHECK_EQUAL(tst_obj.ts_ls_user_emissivity().extent(0), lid_pars.max_surfacewfs);
  BOOST_CHECK_EQUAL(tst_obj.ts_ls_user_emissivity().extent(1), lid_pars.maxstokes);
  BOOST_CHECK_EQUAL(tst_obj.ts_ls_user_emissivity().extent(2), lid_pars.max_user_streams);
  

  // Test initialization
  blitz::Array<double, 5> ts_ls_exactdb_brdfunc_exp(tst_obj.ts_ls_exactdb_brdfunc().shape());
  ts_ls_exactdb_brdfunc_exp = 0;
  BOOST_CHECK_MATRIX_CLOSE_TOL(tst_obj.ts_ls_exactdb_brdfunc(), ts_ls_exactdb_brdfunc_exp, 1e-10);
  blitz::Array<double, 5> ts_ls_brdf_f_0_exp(tst_obj.ts_ls_brdf_f_0().shape());
  ts_ls_brdf_f_0_exp = 0;
  BOOST_CHECK_MATRIX_CLOSE_TOL(tst_obj.ts_ls_brdf_f_0(), ts_ls_brdf_f_0_exp, 1e-10);
  blitz::Array<double, 5> ts_ls_brdf_f_exp(tst_obj.ts_ls_brdf_f().shape());
  ts_ls_brdf_f_exp = 0;
  BOOST_CHECK_MATRIX_CLOSE_TOL(tst_obj.ts_ls_brdf_f(), ts_ls_brdf_f_exp, 1e-10);
  blitz::Array<double, 5> ts_ls_user_brdf_f_0_exp(tst_obj.ts_ls_user_brdf_f_0().shape());
  ts_ls_user_brdf_f_0_exp = 0;
  BOOST_CHECK_MATRIX_CLOSE_TOL(tst_obj.ts_ls_user_brdf_f_0(), ts_ls_user_brdf_f_0_exp, 1e-10);
  blitz::Array<double, 5> ts_ls_user_brdf_f_exp(tst_obj.ts_ls_user_brdf_f().shape());
  ts_ls_user_brdf_f_exp = 0;
  BOOST_CHECK_MATRIX_CLOSE_TOL(tst_obj.ts_ls_user_brdf_f(), ts_ls_user_brdf_f_exp, 1e-10);
  blitz::Array<double, 3> ts_ls_emissivity_exp(tst_obj.ts_ls_emissivity().shape());
  ts_ls_emissivity_exp = 0;
  BOOST_CHECK_MATRIX_CLOSE_TOL(tst_obj.ts_ls_emissivity(), ts_ls_emissivity_exp, 1e-10);
  blitz::Array<double, 3> ts_ls_user_emissivity_exp(tst_obj.ts_ls_user_emissivity().shape());
  ts_ls_user_emissivity_exp = 0;
  BOOST_CHECK_MATRIX_CLOSE_TOL(tst_obj.ts_ls_user_emissivity(), ts_ls_user_emissivity_exp, 1e-10);
  
}

BOOST_AUTO_TEST_CASE(vlidort_linsup_sleave)
{
  // Used for checking dimensions
  VLidort_Pars lid_pars = VLidort_Pars::instance();

  // Test constructor
  VLidort_Linsup_Sleave tst_obj = VLidort_Linsup_Sleave();

  // Test memory space sizes
  tst_obj.check_byte_sizes(); // throws exception on mismatch

  // Test variable shapes
  BOOST_CHECK_EQUAL(tst_obj.ts_lssl_slterm_isotropic().extent(0), lid_pars.max_sleavewfs);
  BOOST_CHECK_EQUAL(tst_obj.ts_lssl_slterm_isotropic().extent(1), lid_pars.maxstokes);
  BOOST_CHECK_EQUAL(tst_obj.ts_lssl_slterm_isotropic().extent(2), lid_pars.maxbeams);
  BOOST_CHECK_EQUAL(tst_obj.ts_lssl_slterm_userangles().extent(0), lid_pars.max_sleavewfs);
  BOOST_CHECK_EQUAL(tst_obj.ts_lssl_slterm_userangles().extent(1), lid_pars.maxstokes);
  BOOST_CHECK_EQUAL(tst_obj.ts_lssl_slterm_userangles().extent(2), lid_pars.max_user_streams);
  BOOST_CHECK_EQUAL(tst_obj.ts_lssl_slterm_userangles().extent(3), lid_pars.max_user_relazms);
  BOOST_CHECK_EQUAL(tst_obj.ts_lssl_slterm_userangles().extent(4), lid_pars.maxbeams);
  BOOST_CHECK_EQUAL(tst_obj.ts_lssl_slterm_f_0().extent(0), lid_pars.max_sleavewfs);
  BOOST_CHECK_EQUAL(tst_obj.ts_lssl_slterm_f_0().extent(1), lid_pars.maxmoments+1);
  BOOST_CHECK_EQUAL(tst_obj.ts_lssl_slterm_f_0().extent(2), lid_pars.maxstokes);
  BOOST_CHECK_EQUAL(tst_obj.ts_lssl_slterm_f_0().extent(3), lid_pars.maxstreams);
  BOOST_CHECK_EQUAL(tst_obj.ts_lssl_slterm_f_0().extent(4), lid_pars.maxbeams);
  BOOST_CHECK_EQUAL(tst_obj.ts_lssl_user_slterm_f_0().extent(0), lid_pars.max_sleavewfs);
  BOOST_CHECK_EQUAL(tst_obj.ts_lssl_user_slterm_f_0().extent(1), lid_pars.maxmoments+1);
  BOOST_CHECK_EQUAL(tst_obj.ts_lssl_user_slterm_f_0().extent(2), lid_pars.maxstokes);
  BOOST_CHECK_EQUAL(tst_obj.ts_lssl_user_slterm_f_0().extent(3), lid_pars.max_user_streams);
  BOOST_CHECK_EQUAL(tst_obj.ts_lssl_user_slterm_f_0().extent(4), lid_pars.maxbeams);
  

  // Test initialization
  blitz::Array<double, 3> ts_lssl_slterm_isotropic_exp(tst_obj.ts_lssl_slterm_isotropic().shape());
  ts_lssl_slterm_isotropic_exp = 0;
  BOOST_CHECK_MATRIX_CLOSE_TOL(tst_obj.ts_lssl_slterm_isotropic(), ts_lssl_slterm_isotropic_exp, 1e-10);
  blitz::Array<double, 5> ts_lssl_slterm_userangles_exp(tst_obj.ts_lssl_slterm_userangles().shape());
  ts_lssl_slterm_userangles_exp = 0;
  BOOST_CHECK_MATRIX_CLOSE_TOL(tst_obj.ts_lssl_slterm_userangles(), ts_lssl_slterm_userangles_exp, 1e-10);
  blitz::Array<double, 5> ts_lssl_slterm_f_0_exp(tst_obj.ts_lssl_slterm_f_0().shape());
  ts_lssl_slterm_f_0_exp = 0;
  BOOST_CHECK_MATRIX_CLOSE_TOL(tst_obj.ts_lssl_slterm_f_0(), ts_lssl_slterm_f_0_exp, 1e-10);
  blitz::Array<double, 5> ts_lssl_user_slterm_f_0_exp(tst_obj.ts_lssl_user_slterm_f_0().shape());
  ts_lssl_user_slterm_f_0_exp = 0;
  BOOST_CHECK_MATRIX_CLOSE_TOL(tst_obj.ts_lssl_user_slterm_f_0(), ts_lssl_user_slterm_f_0_exp, 1e-10);
  
}

BOOST_AUTO_TEST_CASE(vlidort_linsup_ss_col)
{
  // Used for checking dimensions
  VLidort_Pars lid_pars = VLidort_Pars::instance();

  // Test constructor
  VLidort_Linsup_Ss_Col tst_obj = VLidort_Linsup_Ss_Col();

  // Test memory space sizes
  tst_obj.check_byte_sizes(); // throws exception on mismatch

  // Test variable shapes
  BOOST_CHECK_EQUAL(tst_obj.ts_columnwf_ss().extent(0), lid_pars.max_atmoswfs);
  BOOST_CHECK_EQUAL(tst_obj.ts_columnwf_ss().extent(1), lid_pars.max_user_levels);
  BOOST_CHECK_EQUAL(tst_obj.ts_columnwf_ss().extent(2), lid_pars.max_geometries);
  BOOST_CHECK_EQUAL(tst_obj.ts_columnwf_ss().extent(3), lid_pars.maxstokes);
  BOOST_CHECK_EQUAL(tst_obj.ts_columnwf_ss().extent(4), lid_pars.max_directions);
  BOOST_CHECK_EQUAL(tst_obj.ts_columnwf_db().extent(0), lid_pars.max_atmoswfs);
  BOOST_CHECK_EQUAL(tst_obj.ts_columnwf_db().extent(1), lid_pars.max_user_levels);
  BOOST_CHECK_EQUAL(tst_obj.ts_columnwf_db().extent(2), lid_pars.max_geometries);
  BOOST_CHECK_EQUAL(tst_obj.ts_columnwf_db().extent(3), lid_pars.maxstokes);
  

  // Test initialization
  blitz::Array<double, 5> ts_columnwf_ss_exp(tst_obj.ts_columnwf_ss().shape());
  ts_columnwf_ss_exp = 0;
  BOOST_CHECK_MATRIX_CLOSE_TOL(tst_obj.ts_columnwf_ss(), ts_columnwf_ss_exp, 1e-10);
  blitz::Array<double, 4> ts_columnwf_db_exp(tst_obj.ts_columnwf_db().shape());
  ts_columnwf_db_exp = 0;
  BOOST_CHECK_MATRIX_CLOSE_TOL(tst_obj.ts_columnwf_db(), ts_columnwf_db_exp, 1e-10);
  
}

BOOST_AUTO_TEST_CASE(vlidort_linsup_ss_prof)
{
  // Used for checking dimensions
  VLidort_Pars lid_pars = VLidort_Pars::instance();

  // Test constructor
  VLidort_Linsup_Ss_Prof tst_obj = VLidort_Linsup_Ss_Prof();

  // Test memory space sizes
  tst_obj.check_byte_sizes(); // throws exception on mismatch

  // Test variable shapes
  BOOST_CHECK_EQUAL(tst_obj.ts_profilewf_ss().extent(0), lid_pars.max_atmoswfs);
  BOOST_CHECK_EQUAL(tst_obj.ts_profilewf_ss().extent(1), lid_pars.maxlayers);
  BOOST_CHECK_EQUAL(tst_obj.ts_profilewf_ss().extent(2), lid_pars.max_user_levels);
  BOOST_CHECK_EQUAL(tst_obj.ts_profilewf_ss().extent(3), lid_pars.max_geometries);
  BOOST_CHECK_EQUAL(tst_obj.ts_profilewf_ss().extent(4), lid_pars.maxstokes);
  BOOST_CHECK_EQUAL(tst_obj.ts_profilewf_ss().extent(5), lid_pars.max_directions);
  BOOST_CHECK_EQUAL(tst_obj.ts_profilewf_db().extent(0), lid_pars.max_atmoswfs);
  BOOST_CHECK_EQUAL(tst_obj.ts_profilewf_db().extent(1), lid_pars.maxlayers);
  BOOST_CHECK_EQUAL(tst_obj.ts_profilewf_db().extent(2), lid_pars.max_user_levels);
  BOOST_CHECK_EQUAL(tst_obj.ts_profilewf_db().extent(3), lid_pars.max_geometries);
  BOOST_CHECK_EQUAL(tst_obj.ts_profilewf_db().extent(4), lid_pars.maxstokes);
  

  // Test initialization
  blitz::Array<double, 6> ts_profilewf_ss_exp(tst_obj.ts_profilewf_ss().shape());
  ts_profilewf_ss_exp = 0;
  BOOST_CHECK_MATRIX_CLOSE_TOL(tst_obj.ts_profilewf_ss(), ts_profilewf_ss_exp, 1e-10);
  blitz::Array<double, 5> ts_profilewf_db_exp(tst_obj.ts_profilewf_db().shape());
  ts_profilewf_db_exp = 0;
  BOOST_CHECK_MATRIX_CLOSE_TOL(tst_obj.ts_profilewf_db(), ts_profilewf_db_exp, 1e-10);
  
}

BOOST_AUTO_TEST_CASE(vlidort_linsup_ss_surf)
{
  // Used for checking dimensions
  VLidort_Pars lid_pars = VLidort_Pars::instance();

  // Test constructor
  VLidort_Linsup_Ss_Surf tst_obj = VLidort_Linsup_Ss_Surf();

  // Test memory space sizes
  tst_obj.check_byte_sizes(); // throws exception on mismatch

  // Test variable shapes
  BOOST_CHECK_EQUAL(tst_obj.ts_surfacewf_db().extent(0), lid_pars.max_surfacewfs);
  BOOST_CHECK_EQUAL(tst_obj.ts_surfacewf_db().extent(1), lid_pars.max_user_levels);
  BOOST_CHECK_EQUAL(tst_obj.ts_surfacewf_db().extent(2), lid_pars.max_geometries);
  BOOST_CHECK_EQUAL(tst_obj.ts_surfacewf_db().extent(3), lid_pars.maxstokes);
  

  // Test initialization
  blitz::Array<double, 4> ts_surfacewf_db_exp(tst_obj.ts_surfacewf_db().shape());
  ts_surfacewf_db_exp = 0;
  BOOST_CHECK_MATRIX_CLOSE_TOL(tst_obj.ts_surfacewf_db(), ts_surfacewf_db_exp, 1e-10);
  
}

BOOST_AUTO_TEST_CASE(vlidort_linsup_ss)
{
  // Used for checking dimensions
  VLidort_Pars lid_pars = VLidort_Pars::instance();

  // Test constructor
  VLidort_Linsup_Ss tst_obj = VLidort_Linsup_Ss();

  // Test memory space sizes
  tst_obj.check_byte_sizes(); // throws exception on mismatch

  // Test variable shapes
  

  // Test initialization
  
}

BOOST_AUTO_TEST_CASE(vlidort_linsup_inout)
{
  // Used for checking dimensions
  VLidort_Pars lid_pars = VLidort_Pars::instance();

  // Test constructor
  VLidort_Linsup_Inout tst_obj = VLidort_Linsup_Inout();

  // Test memory space sizes
  tst_obj.check_byte_sizes(); // throws exception on mismatch

  // Test variable shapes
  

  // Test initialization
  
}

BOOST_AUTO_TEST_CASE(vlidort_main_outputs)
{
  // Used for checking dimensions
  VLidort_Pars lid_pars = VLidort_Pars::instance();

  // Test constructor
  VLidort_Main_Outputs tst_obj = VLidort_Main_Outputs();

  // Test memory space sizes
  tst_obj.check_byte_sizes(); // throws exception on mismatch

  // Test variable shapes
  BOOST_CHECK_EQUAL(tst_obj.ts_stokes().extent(0), lid_pars.max_user_levels);
  BOOST_CHECK_EQUAL(tst_obj.ts_stokes().extent(1), lid_pars.max_geometries);
  BOOST_CHECK_EQUAL(tst_obj.ts_stokes().extent(2), lid_pars.maxstokes);
  BOOST_CHECK_EQUAL(tst_obj.ts_stokes().extent(3), lid_pars.max_directions);
  BOOST_CHECK_EQUAL(tst_obj.ts_meanst_diffuse().extent(0), lid_pars.max_user_levels);
  BOOST_CHECK_EQUAL(tst_obj.ts_meanst_diffuse().extent(1), lid_pars.max_szangles);
  BOOST_CHECK_EQUAL(tst_obj.ts_meanst_diffuse().extent(2), lid_pars.maxstokes);
  BOOST_CHECK_EQUAL(tst_obj.ts_meanst_diffuse().extent(3), lid_pars.max_directions);
  BOOST_CHECK_EQUAL(tst_obj.ts_flux_diffuse().extent(0), lid_pars.max_user_levels);
  BOOST_CHECK_EQUAL(tst_obj.ts_flux_diffuse().extent(1), lid_pars.max_szangles);
  BOOST_CHECK_EQUAL(tst_obj.ts_flux_diffuse().extent(2), lid_pars.maxstokes);
  BOOST_CHECK_EQUAL(tst_obj.ts_flux_diffuse().extent(3), lid_pars.max_directions);
  BOOST_CHECK_EQUAL(tst_obj.ts_dnmeanst_direct().extent(0), lid_pars.max_user_levels);
  BOOST_CHECK_EQUAL(tst_obj.ts_dnmeanst_direct().extent(1), lid_pars.max_szangles);
  BOOST_CHECK_EQUAL(tst_obj.ts_dnmeanst_direct().extent(2), lid_pars.maxstokes);
  BOOST_CHECK_EQUAL(tst_obj.ts_dnflux_direct().extent(0), lid_pars.max_user_levels);
  BOOST_CHECK_EQUAL(tst_obj.ts_dnflux_direct().extent(1), lid_pars.max_szangles);
  BOOST_CHECK_EQUAL(tst_obj.ts_dnflux_direct().extent(2), lid_pars.maxstokes);
  BOOST_CHECK_EQUAL(tst_obj.ts_albmed_user().extent(0), lid_pars.maxstokes);
  BOOST_CHECK_EQUAL(tst_obj.ts_albmed_user().extent(1), lid_pars.max_user_vzangles);
  BOOST_CHECK_EQUAL(tst_obj.ts_trnmed_user().extent(0), lid_pars.maxstokes);
  BOOST_CHECK_EQUAL(tst_obj.ts_trnmed_user().extent(1), lid_pars.max_user_vzangles);
  BOOST_CHECK_EQUAL(tst_obj.ts_albmed_fluxes().extent(0), lid_pars.maxstokes);
  BOOST_CHECK_EQUAL(tst_obj.ts_albmed_fluxes().extent(1), lid_pars.2);
  BOOST_CHECK_EQUAL(tst_obj.ts_trnmed_fluxes().extent(0), lid_pars.maxstokes);
  BOOST_CHECK_EQUAL(tst_obj.ts_trnmed_fluxes().extent(1), lid_pars.2);
  BOOST_CHECK_EQUAL(tst_obj.ts_planetary_transterm().extent(0), lid_pars.maxstokes);
  BOOST_CHECK_EQUAL(tst_obj.ts_planetary_transterm().extent(1), lid_pars.max_geometries);
  BOOST_CHECK_EQUAL(tst_obj.ts_pathgeoms().extent(0), lid_pars.2);
  BOOST_CHECK_EQUAL(tst_obj.ts_pathgeoms().extent(1), lid_pars.maxlayers+1);
  BOOST_CHECK_EQUAL(tst_obj.ts_lostrans().extent(0), lid_pars.max_szangles);
  BOOST_CHECK_EQUAL(tst_obj.ts_lostrans().extent(1), lid_pars.maxlayers);
  BOOST_CHECK_EQUAL(tst_obj.ts_layer_mssts().extent(0), lid_pars.max_szangles);
  BOOST_CHECK_EQUAL(tst_obj.ts_layer_mssts().extent(1), lid_pars.maxstokes);
  BOOST_CHECK_EQUAL(tst_obj.ts_layer_mssts().extent(2), lid_pars.maxlayers);
  BOOST_CHECK_EQUAL(tst_obj.ts_surf_mssts().extent(0), lid_pars.max_szangles);
  BOOST_CHECK_EQUAL(tst_obj.ts_surf_mssts().extent(1), lid_pars.maxstokes);
  BOOST_CHECK_EQUAL(tst_obj.ts_contribs().extent(0), lid_pars.max_geometries);
  BOOST_CHECK_EQUAL(tst_obj.ts_contribs().extent(1), lid_pars.maxstokes);
  BOOST_CHECK_EQUAL(tst_obj.ts_contribs().extent(2), lid_pars.maxlayers);
  BOOST_CHECK_EQUAL(tst_obj.ts_fourier_saved().extent(0), lid_pars.max_szangles);
  BOOST_CHECK_EQUAL(tst_obj.ts_sza_offsets().extent(0), lid_pars.max_szangles);
  BOOST_CHECK_EQUAL(tst_obj.ts_vza_offsets().extent(0), lid_pars.max_szangles);
  BOOST_CHECK_EQUAL(tst_obj.ts_vza_offsets().extent(1), lid_pars.max_user_vzangles);
  BOOST_CHECK_EQUAL(tst_obj.ts_szd_offsets().extent(0), lid_pars.max_szangles);
  BOOST_CHECK_EQUAL(tst_obj.ts_solarbeam_boatrans().extent(0), lid_pars.maxbeams);
  

  // Test initialization
  blitz::Array<double, 4> ts_stokes_exp(tst_obj.ts_stokes().shape());
  ts_stokes_exp = 0;
  BOOST_CHECK_MATRIX_CLOSE_TOL(tst_obj.ts_stokes(), ts_stokes_exp, 1e-10);
  blitz::Array<double, 4> ts_meanst_diffuse_exp(tst_obj.ts_meanst_diffuse().shape());
  ts_meanst_diffuse_exp = 0;
  BOOST_CHECK_MATRIX_CLOSE_TOL(tst_obj.ts_meanst_diffuse(), ts_meanst_diffuse_exp, 1e-10);
  blitz::Array<double, 4> ts_flux_diffuse_exp(tst_obj.ts_flux_diffuse().shape());
  ts_flux_diffuse_exp = 0;
  BOOST_CHECK_MATRIX_CLOSE_TOL(tst_obj.ts_flux_diffuse(), ts_flux_diffuse_exp, 1e-10);
  blitz::Array<double, 3> ts_dnmeanst_direct_exp(tst_obj.ts_dnmeanst_direct().shape());
  ts_dnmeanst_direct_exp = 0;
  BOOST_CHECK_MATRIX_CLOSE_TOL(tst_obj.ts_dnmeanst_direct(), ts_dnmeanst_direct_exp, 1e-10);
  blitz::Array<double, 3> ts_dnflux_direct_exp(tst_obj.ts_dnflux_direct().shape());
  ts_dnflux_direct_exp = 0;
  BOOST_CHECK_MATRIX_CLOSE_TOL(tst_obj.ts_dnflux_direct(), ts_dnflux_direct_exp, 1e-10);
  blitz::Array<double, 2> ts_albmed_user_exp(tst_obj.ts_albmed_user().shape());
  ts_albmed_user_exp = 0;
  BOOST_CHECK_MATRIX_CLOSE_TOL(tst_obj.ts_albmed_user(), ts_albmed_user_exp, 1e-10);
  blitz::Array<double, 2> ts_trnmed_user_exp(tst_obj.ts_trnmed_user().shape());
  ts_trnmed_user_exp = 0;
  BOOST_CHECK_MATRIX_CLOSE_TOL(tst_obj.ts_trnmed_user(), ts_trnmed_user_exp, 1e-10);
  blitz::Array<double, 2> ts_albmed_fluxes_exp(tst_obj.ts_albmed_fluxes().shape());
  ts_albmed_fluxes_exp = 0;
  BOOST_CHECK_MATRIX_CLOSE_TOL(tst_obj.ts_albmed_fluxes(), ts_albmed_fluxes_exp, 1e-10);
  blitz::Array<double, 2> ts_trnmed_fluxes_exp(tst_obj.ts_trnmed_fluxes().shape());
  ts_trnmed_fluxes_exp = 0;
  BOOST_CHECK_MATRIX_CLOSE_TOL(tst_obj.ts_trnmed_fluxes(), ts_trnmed_fluxes_exp, 1e-10);
  blitz::Array<double, 2> ts_planetary_transterm_exp(tst_obj.ts_planetary_transterm().shape());
  ts_planetary_transterm_exp = 0;
  BOOST_CHECK_MATRIX_CLOSE_TOL(tst_obj.ts_planetary_transterm(), ts_planetary_transterm_exp, 1e-10);
  BOOST_CHECK_CLOSE(tst_obj.ts_planetary_sbterm(), 0, 1e-10);
  blitz::Array<double, 2> ts_pathgeoms_exp(tst_obj.ts_pathgeoms().shape());
  ts_pathgeoms_exp = 0;
  BOOST_CHECK_MATRIX_CLOSE_TOL(tst_obj.ts_pathgeoms(), ts_pathgeoms_exp, 1e-10);
  blitz::Array<double, 2> ts_lostrans_exp(tst_obj.ts_lostrans().shape());
  ts_lostrans_exp = 0;
  BOOST_CHECK_MATRIX_CLOSE_TOL(tst_obj.ts_lostrans(), ts_lostrans_exp, 1e-10);
  blitz::Array<double, 3> ts_layer_mssts_exp(tst_obj.ts_layer_mssts().shape());
  ts_layer_mssts_exp = 0;
  BOOST_CHECK_MATRIX_CLOSE_TOL(tst_obj.ts_layer_mssts(), ts_layer_mssts_exp, 1e-10);
  blitz::Array<double, 2> ts_surf_mssts_exp(tst_obj.ts_surf_mssts().shape());
  ts_surf_mssts_exp = 0;
  BOOST_CHECK_MATRIX_CLOSE_TOL(tst_obj.ts_surf_mssts(), ts_surf_mssts_exp, 1e-10);
  blitz::Array<double, 3> ts_contribs_exp(tst_obj.ts_contribs().shape());
  ts_contribs_exp = 0;
  BOOST_CHECK_MATRIX_CLOSE_TOL(tst_obj.ts_contribs(), ts_contribs_exp, 1e-10);
  blitz::Array<int, 1> ts_fourier_saved_exp(tst_obj.ts_fourier_saved().shape());
  ts_fourier_saved_exp = 0;
  BOOST_CHECK_MATRIX_CLOSE_TOL(tst_obj.ts_fourier_saved(), ts_fourier_saved_exp, 1e-10);
  BOOST_CHECK_EQUAL(tst_obj.ts_n_geometries(), 0);
  blitz::Array<int, 1> ts_sza_offsets_exp(tst_obj.ts_sza_offsets().shape());
  ts_sza_offsets_exp = 0;
  BOOST_CHECK_MATRIX_CLOSE_TOL(tst_obj.ts_sza_offsets(), ts_sza_offsets_exp, 1e-10);
  blitz::Array<int, 2> ts_vza_offsets_exp(tst_obj.ts_vza_offsets().shape());
  ts_vza_offsets_exp = 0;
  BOOST_CHECK_MATRIX_CLOSE_TOL(tst_obj.ts_vza_offsets(), ts_vza_offsets_exp, 1e-10);
  blitz::Array<int, 1> ts_szd_offsets_exp(tst_obj.ts_szd_offsets().shape());
  ts_szd_offsets_exp = 0;
  BOOST_CHECK_MATRIX_CLOSE_TOL(tst_obj.ts_szd_offsets(), ts_szd_offsets_exp, 1e-10);
  blitz::Array<double, 1> ts_solarbeam_boatrans_exp(tst_obj.ts_solarbeam_boatrans().shape());
  ts_solarbeam_boatrans_exp = 0;
  BOOST_CHECK_MATRIX_CLOSE_TOL(tst_obj.ts_solarbeam_boatrans(), ts_solarbeam_boatrans_exp, 1e-10);
  
}

BOOST_AUTO_TEST_CASE(vlidort_wladjusted_outputs)
{
  // Used for checking dimensions
  VLidort_Pars lid_pars = VLidort_Pars::instance();

  // Test constructor
  VLidort_Wladjusted_Outputs tst_obj = VLidort_Wladjusted_Outputs();

  // Test memory space sizes
  tst_obj.check_byte_sizes(); // throws exception on mismatch

  // Test variable shapes
  BOOST_CHECK_EQUAL(tst_obj.ts_wladjusted_isotropic().extent(0), lid_pars.maxstokes);
  BOOST_CHECK_EQUAL(tst_obj.ts_wladjusted_isotropic().extent(1), lid_pars.maxbeams);
  BOOST_CHECK_EQUAL(tst_obj.ts_wladjusted_direct().extent(0), lid_pars.maxstokes);
  BOOST_CHECK_EQUAL(tst_obj.ts_wladjusted_direct().extent(1), lid_pars.max_user_vzangles);
  BOOST_CHECK_EQUAL(tst_obj.ts_wladjusted_direct().extent(2), lid_pars.max_user_relazms);
  BOOST_CHECK_EQUAL(tst_obj.ts_wladjusted_direct().extent(3), lid_pars.maxbeams);
  BOOST_CHECK_EQUAL(tst_obj.ts_wladjusted_f_ords_0().extent(0), lid_pars.maxmoments+1);
  BOOST_CHECK_EQUAL(tst_obj.ts_wladjusted_f_ords_0().extent(1), lid_pars.maxstokes);
  BOOST_CHECK_EQUAL(tst_obj.ts_wladjusted_f_ords_0().extent(2), lid_pars.maxstreams);
  BOOST_CHECK_EQUAL(tst_obj.ts_wladjusted_f_ords_0().extent(3), lid_pars.maxbeams);
  BOOST_CHECK_EQUAL(tst_obj.ts_wladjusted_f_user_0().extent(0), lid_pars.maxmoments+1);
  BOOST_CHECK_EQUAL(tst_obj.ts_wladjusted_f_user_0().extent(1), lid_pars.maxstokes);
  BOOST_CHECK_EQUAL(tst_obj.ts_wladjusted_f_user_0().extent(2), lid_pars.max_user_vzangles);
  BOOST_CHECK_EQUAL(tst_obj.ts_wladjusted_f_user_0().extent(3), lid_pars.maxbeams);
  

  // Test initialization
  blitz::Array<double, 2> ts_wladjusted_isotropic_exp(tst_obj.ts_wladjusted_isotropic().shape());
  ts_wladjusted_isotropic_exp = 0;
  BOOST_CHECK_MATRIX_CLOSE_TOL(tst_obj.ts_wladjusted_isotropic(), ts_wladjusted_isotropic_exp, 1e-10);
  blitz::Array<double, 4> ts_wladjusted_direct_exp(tst_obj.ts_wladjusted_direct().shape());
  ts_wladjusted_direct_exp = 0;
  BOOST_CHECK_MATRIX_CLOSE_TOL(tst_obj.ts_wladjusted_direct(), ts_wladjusted_direct_exp, 1e-10);
  blitz::Array<double, 4> ts_wladjusted_f_ords_0_exp(tst_obj.ts_wladjusted_f_ords_0().shape());
  ts_wladjusted_f_ords_0_exp = 0;
  BOOST_CHECK_MATRIX_CLOSE_TOL(tst_obj.ts_wladjusted_f_ords_0(), ts_wladjusted_f_ords_0_exp, 1e-10);
  blitz::Array<double, 4> ts_wladjusted_f_user_0_exp(tst_obj.ts_wladjusted_f_user_0().shape());
  ts_wladjusted_f_user_0_exp = 0;
  BOOST_CHECK_MATRIX_CLOSE_TOL(tst_obj.ts_wladjusted_f_user_0(), ts_wladjusted_f_user_0_exp, 1e-10);
  
}

BOOST_AUTO_TEST_CASE(vlidort_exception_handling)
{
  // Used for checking dimensions
  VLidort_Pars lid_pars = VLidort_Pars::instance();

  // Test constructor
  VLidort_Exception_Handling tst_obj = VLidort_Exception_Handling();

  // Test memory space sizes
  tst_obj.check_byte_sizes(); // throws exception on mismatch

  // Test variable shapes
  

  // Test initialization
  BOOST_CHECK_EQUAL(tst_obj.ts_status_inputcheck(), 0);
  BOOST_CHECK_EQUAL(tst_obj.ts_ncheckmessages(), 0);
  BOOST_CHECK_EQUAL(tst_obj.ts_status_calculation(), 0);
  
}

BOOST_AUTO_TEST_CASE(vlidort_input_exception_handling)
{
  // Used for checking dimensions
  VLidort_Pars lid_pars = VLidort_Pars::instance();

  // Test constructor
  VLidort_Input_Exception_Handling tst_obj = VLidort_Input_Exception_Handling();

  // Test memory space sizes
  tst_obj.check_byte_sizes(); // throws exception on mismatch

  // Test variable shapes
  

  // Test initialization
  BOOST_CHECK_EQUAL(tst_obj.ts_status_inputread(), 0);
  BOOST_CHECK_EQUAL(tst_obj.ts_ninputmessages(), 0);
  
}

BOOST_AUTO_TEST_CASE(vlidort_outputs)
{
  // Used for checking dimensions
  VLidort_Pars lid_pars = VLidort_Pars::instance();

  // Test constructor
  VLidort_Outputs tst_obj = VLidort_Outputs();

  // Test memory space sizes
  tst_obj.check_byte_sizes(); // throws exception on mismatch

  // Test variable shapes
  

  // Test initialization
  
}

BOOST_AUTO_TEST_CASE(vlidort_sup_brdf)
{
  // Used for checking dimensions
  VLidort_Pars lid_pars = VLidort_Pars::instance();

  // Test constructor
  VLidort_Sup_Brdf tst_obj = VLidort_Sup_Brdf();

  // Test memory space sizes
  tst_obj.check_byte_sizes(); // throws exception on mismatch

  // Test variable shapes
  BOOST_CHECK_EQUAL(tst_obj.ts_exactdb_brdfunc().extent(0), lid_pars.maxstokes_sq);
  BOOST_CHECK_EQUAL(tst_obj.ts_exactdb_brdfunc().extent(1), lid_pars.max_user_streams);
  BOOST_CHECK_EQUAL(tst_obj.ts_exactdb_brdfunc().extent(2), lid_pars.max_user_relazms);
  BOOST_CHECK_EQUAL(tst_obj.ts_exactdb_brdfunc().extent(3), lid_pars.maxbeams);
  BOOST_CHECK_EQUAL(tst_obj.ts_brdf_f_0().extent(0), lid_pars.maxmoments+1);
  BOOST_CHECK_EQUAL(tst_obj.ts_brdf_f_0().extent(1), lid_pars.maxstokes_sq);
  BOOST_CHECK_EQUAL(tst_obj.ts_brdf_f_0().extent(2), lid_pars.maxstreams);
  BOOST_CHECK_EQUAL(tst_obj.ts_brdf_f_0().extent(3), lid_pars.maxbeams);
  BOOST_CHECK_EQUAL(tst_obj.ts_brdf_f().extent(0), lid_pars.maxmoments+1);
  BOOST_CHECK_EQUAL(tst_obj.ts_brdf_f().extent(1), lid_pars.maxstokes_sq);
  BOOST_CHECK_EQUAL(tst_obj.ts_brdf_f().extent(2), lid_pars.maxstreams);
  BOOST_CHECK_EQUAL(tst_obj.ts_brdf_f().extent(3), lid_pars.maxstreams);
  BOOST_CHECK_EQUAL(tst_obj.ts_user_brdf_f_0().extent(0), lid_pars.maxmoments+1);
  BOOST_CHECK_EQUAL(tst_obj.ts_user_brdf_f_0().extent(1), lid_pars.maxstokes_sq);
  BOOST_CHECK_EQUAL(tst_obj.ts_user_brdf_f_0().extent(2), lid_pars.max_user_streams);
  BOOST_CHECK_EQUAL(tst_obj.ts_user_brdf_f_0().extent(3), lid_pars.maxbeams);
  BOOST_CHECK_EQUAL(tst_obj.ts_user_brdf_f().extent(0), lid_pars.maxmoments+1);
  BOOST_CHECK_EQUAL(tst_obj.ts_user_brdf_f().extent(1), lid_pars.maxstokes_sq);
  BOOST_CHECK_EQUAL(tst_obj.ts_user_brdf_f().extent(2), lid_pars.max_user_streams);
  BOOST_CHECK_EQUAL(tst_obj.ts_user_brdf_f().extent(3), lid_pars.maxstreams);
  BOOST_CHECK_EQUAL(tst_obj.ts_emissivity().extent(0), lid_pars.maxstokes);
  BOOST_CHECK_EQUAL(tst_obj.ts_emissivity().extent(1), lid_pars.maxstreams);
  BOOST_CHECK_EQUAL(tst_obj.ts_user_emissivity().extent(0), lid_pars.maxstokes);
  BOOST_CHECK_EQUAL(tst_obj.ts_user_emissivity().extent(1), lid_pars.max_user_streams);
  

  // Test initialization
  blitz::Array<double, 4> ts_exactdb_brdfunc_exp(tst_obj.ts_exactdb_brdfunc().shape());
  ts_exactdb_brdfunc_exp = 0;
  BOOST_CHECK_MATRIX_CLOSE_TOL(tst_obj.ts_exactdb_brdfunc(), ts_exactdb_brdfunc_exp, 1e-10);
  blitz::Array<double, 4> ts_brdf_f_0_exp(tst_obj.ts_brdf_f_0().shape());
  ts_brdf_f_0_exp = 0;
  BOOST_CHECK_MATRIX_CLOSE_TOL(tst_obj.ts_brdf_f_0(), ts_brdf_f_0_exp, 1e-10);
  blitz::Array<double, 4> ts_brdf_f_exp(tst_obj.ts_brdf_f().shape());
  ts_brdf_f_exp = 0;
  BOOST_CHECK_MATRIX_CLOSE_TOL(tst_obj.ts_brdf_f(), ts_brdf_f_exp, 1e-10);
  blitz::Array<double, 4> ts_user_brdf_f_0_exp(tst_obj.ts_user_brdf_f_0().shape());
  ts_user_brdf_f_0_exp = 0;
  BOOST_CHECK_MATRIX_CLOSE_TOL(tst_obj.ts_user_brdf_f_0(), ts_user_brdf_f_0_exp, 1e-10);
  blitz::Array<double, 4> ts_user_brdf_f_exp(tst_obj.ts_user_brdf_f().shape());
  ts_user_brdf_f_exp = 0;
  BOOST_CHECK_MATRIX_CLOSE_TOL(tst_obj.ts_user_brdf_f(), ts_user_brdf_f_exp, 1e-10);
  blitz::Array<double, 2> ts_emissivity_exp(tst_obj.ts_emissivity().shape());
  ts_emissivity_exp = 0;
  BOOST_CHECK_MATRIX_CLOSE_TOL(tst_obj.ts_emissivity(), ts_emissivity_exp, 1e-10);
  blitz::Array<double, 2> ts_user_emissivity_exp(tst_obj.ts_user_emissivity().shape());
  ts_user_emissivity_exp = 0;
  BOOST_CHECK_MATRIX_CLOSE_TOL(tst_obj.ts_user_emissivity(), ts_user_emissivity_exp, 1e-10);
  
}

BOOST_AUTO_TEST_CASE(vlidort_sup_sleave)
{
  // Used for checking dimensions
  VLidort_Pars lid_pars = VLidort_Pars::instance();

  // Test constructor
  VLidort_Sup_Sleave tst_obj = VLidort_Sup_Sleave();

  // Test memory space sizes
  tst_obj.check_byte_sizes(); // throws exception on mismatch

  // Test variable shapes
  BOOST_CHECK_EQUAL(tst_obj.ts_slterm_isotropic().extent(0), lid_pars.maxstokes);
  BOOST_CHECK_EQUAL(tst_obj.ts_slterm_isotropic().extent(1), lid_pars.maxbeams);
  BOOST_CHECK_EQUAL(tst_obj.ts_slterm_userangles().extent(0), lid_pars.maxstokes);
  BOOST_CHECK_EQUAL(tst_obj.ts_slterm_userangles().extent(1), lid_pars.max_user_streams);
  BOOST_CHECK_EQUAL(tst_obj.ts_slterm_userangles().extent(2), lid_pars.max_user_relazms);
  BOOST_CHECK_EQUAL(tst_obj.ts_slterm_userangles().extent(3), lid_pars.maxbeams);
  BOOST_CHECK_EQUAL(tst_obj.ts_slterm_f_0().extent(0), lid_pars.maxmoments+1);
  BOOST_CHECK_EQUAL(tst_obj.ts_slterm_f_0().extent(1), lid_pars.maxstokes);
  BOOST_CHECK_EQUAL(tst_obj.ts_slterm_f_0().extent(2), lid_pars.maxstreams);
  BOOST_CHECK_EQUAL(tst_obj.ts_slterm_f_0().extent(3), lid_pars.maxbeams);
  BOOST_CHECK_EQUAL(tst_obj.ts_user_slterm_f_0().extent(0), lid_pars.maxmoments+1);
  BOOST_CHECK_EQUAL(tst_obj.ts_user_slterm_f_0().extent(1), lid_pars.maxstokes);
  BOOST_CHECK_EQUAL(tst_obj.ts_user_slterm_f_0().extent(2), lid_pars.max_user_streams);
  BOOST_CHECK_EQUAL(tst_obj.ts_user_slterm_f_0().extent(3), lid_pars.maxbeams);
  

  // Test initialization
  blitz::Array<double, 2> ts_slterm_isotropic_exp(tst_obj.ts_slterm_isotropic().shape());
  ts_slterm_isotropic_exp = 0;
  BOOST_CHECK_MATRIX_CLOSE_TOL(tst_obj.ts_slterm_isotropic(), ts_slterm_isotropic_exp, 1e-10);
  blitz::Array<double, 4> ts_slterm_userangles_exp(tst_obj.ts_slterm_userangles().shape());
  ts_slterm_userangles_exp = 0;
  BOOST_CHECK_MATRIX_CLOSE_TOL(tst_obj.ts_slterm_userangles(), ts_slterm_userangles_exp, 1e-10);
  blitz::Array<double, 4> ts_slterm_f_0_exp(tst_obj.ts_slterm_f_0().shape());
  ts_slterm_f_0_exp = 0;
  BOOST_CHECK_MATRIX_CLOSE_TOL(tst_obj.ts_slterm_f_0(), ts_slterm_f_0_exp, 1e-10);
  blitz::Array<double, 4> ts_user_slterm_f_0_exp(tst_obj.ts_user_slterm_f_0().shape());
  ts_user_slterm_f_0_exp = 0;
  BOOST_CHECK_MATRIX_CLOSE_TOL(tst_obj.ts_user_slterm_f_0(), ts_user_slterm_f_0_exp, 1e-10);
  
}

BOOST_AUTO_TEST_CASE(vlidort_sup_ss)
{
  // Used for checking dimensions
  VLidort_Pars lid_pars = VLidort_Pars::instance();

  // Test constructor
  VLidort_Sup_Ss tst_obj = VLidort_Sup_Ss();

  // Test memory space sizes
  tst_obj.check_byte_sizes(); // throws exception on mismatch

  // Test variable shapes
  BOOST_CHECK_EQUAL(tst_obj.ts_stokes_ss().extent(0), lid_pars.max_user_levels);
  BOOST_CHECK_EQUAL(tst_obj.ts_stokes_ss().extent(1), lid_pars.max_geometries);
  BOOST_CHECK_EQUAL(tst_obj.ts_stokes_ss().extent(2), lid_pars.maxstokes);
  BOOST_CHECK_EQUAL(tst_obj.ts_stokes_ss().extent(3), lid_pars.max_directions);
  BOOST_CHECK_EQUAL(tst_obj.ts_stokes_db().extent(0), lid_pars.max_user_levels);
  BOOST_CHECK_EQUAL(tst_obj.ts_stokes_db().extent(1), lid_pars.max_geometries);
  BOOST_CHECK_EQUAL(tst_obj.ts_stokes_db().extent(2), lid_pars.maxstokes);
  BOOST_CHECK_EQUAL(tst_obj.ts_contribs_ss().extent(0), lid_pars.max_geometries);
  BOOST_CHECK_EQUAL(tst_obj.ts_contribs_ss().extent(1), lid_pars.maxstokes);
  BOOST_CHECK_EQUAL(tst_obj.ts_contribs_ss().extent(2), lid_pars.maxlayers);
  

  // Test initialization
  blitz::Array<double, 4> ts_stokes_ss_exp(tst_obj.ts_stokes_ss().shape());
  ts_stokes_ss_exp = 0;
  BOOST_CHECK_MATRIX_CLOSE_TOL(tst_obj.ts_stokes_ss(), ts_stokes_ss_exp, 1e-10);
  blitz::Array<double, 3> ts_stokes_db_exp(tst_obj.ts_stokes_db().shape());
  ts_stokes_db_exp = 0;
  BOOST_CHECK_MATRIX_CLOSE_TOL(tst_obj.ts_stokes_db(), ts_stokes_db_exp, 1e-10);
  blitz::Array<double, 3> ts_contribs_ss_exp(tst_obj.ts_contribs_ss().shape());
  ts_contribs_ss_exp = 0;
  BOOST_CHECK_MATRIX_CLOSE_TOL(tst_obj.ts_contribs_ss(), ts_contribs_ss_exp, 1e-10);
  
}

BOOST_AUTO_TEST_CASE(vlidort_sup_inout)
{
  // Used for checking dimensions
  VLidort_Pars lid_pars = VLidort_Pars::instance();

  // Test constructor
  VLidort_Sup_Inout tst_obj = VLidort_Sup_Inout();

  // Test memory space sizes
  tst_obj.check_byte_sizes(); // throws exception on mismatch

  // Test variable shapes
  

  // Test initialization
  
}

BOOST_AUTO_TEST_CASE(vlidort_fixed_boolean)
{
  // Used for checking dimensions
  VLidort_Pars lid_pars = VLidort_Pars::instance();

  // Test constructor
  VLidort_Fixed_Boolean tst_obj = VLidort_Fixed_Boolean();

  // Test memory space sizes
  tst_obj.check_byte_sizes(); // throws exception on mismatch

  // Test variable shapes
  BOOST_CHECK_EQUAL(tst_obj.ts_do_albtrn_media().extent(0), lid_pars.2);
  

  // Test initialization
  BOOST_CHECK_EQUAL(tst_obj.ts_do_fullrad_mode(), false);
  BOOST_CHECK_EQUAL(tst_obj.ts_do_thermal_emission(), false);
  BOOST_CHECK_EQUAL(tst_obj.ts_do_surface_emission(), false);
  BOOST_CHECK_EQUAL(tst_obj.ts_do_plane_parallel(), false);
  BOOST_CHECK_EQUAL(tst_obj.ts_do_upwelling(), false);
  BOOST_CHECK_EQUAL(tst_obj.ts_do_dnwelling(), false);
  BOOST_CHECK_EQUAL(tst_obj.ts_do_toa_contribs(), false);
  BOOST_CHECK_EQUAL(tst_obj.ts_do_lambertian_surface(), false);
  BOOST_CHECK_EQUAL(tst_obj.ts_do_specialist_option_1(), false);
  BOOST_CHECK_EQUAL(tst_obj.ts_do_specialist_option_2(), false);
  BOOST_CHECK_EQUAL(tst_obj.ts_do_specialist_option_3(), false);
  BOOST_CHECK_EQUAL(tst_obj.ts_do_surface_leaving(), false);
  BOOST_CHECK_EQUAL(tst_obj.ts_do_sl_isotropic(), false);
  BOOST_CHECK_EQUAL(tst_obj.ts_do_water_leaving(), false);
  BOOST_CHECK_EQUAL(tst_obj.ts_do_fluorescence(), false);
  BOOST_CHECK_EQUAL(tst_obj.ts_do_tf_iteration(), false);
  BOOST_CHECK_EQUAL(tst_obj.ts_do_wladjusted_output(), false);
  BOOST_CHECK_EQUAL(tst_obj.ts_do_toa_illumination(), false);
  BOOST_CHECK_EQUAL(tst_obj.ts_do_boa_illumination(), false);
  blitz::Array<bool, 1> ts_do_albtrn_media_exp(tst_obj.ts_do_albtrn_media().shape());
  ts_do_albtrn_media_exp = false;
  BOOST_CHECK_MATRIX_CLOSE_TOL(tst_obj.ts_do_albtrn_media(), ts_do_albtrn_media_exp, 1e-10);
  BOOST_CHECK_EQUAL(tst_obj.ts_do_planetary_problem(), false);
  BOOST_CHECK_EQUAL(tst_obj.ts_do_mssts(), false);
  BOOST_CHECK_EQUAL(tst_obj.ts_do_fourier0_nstokes2(), false);
  
}

BOOST_AUTO_TEST_CASE(vlidort_fixed_control)
{
  // Used for checking dimensions
  VLidort_Pars lid_pars = VLidort_Pars::instance();

  // Test constructor
  VLidort_Fixed_Control tst_obj = VLidort_Fixed_Control();

  // Test memory space sizes
  tst_obj.check_byte_sizes(); // throws exception on mismatch

  // Test variable shapes
  

  // Test initialization
  BOOST_CHECK_EQUAL(tst_obj.ts_taylor_order(), 0);
  BOOST_CHECK_EQUAL(tst_obj.ts_nstokes(), 0);
  BOOST_CHECK_EQUAL(tst_obj.ts_nstreams(), 0);
  BOOST_CHECK_EQUAL(tst_obj.ts_nlayers(), 0);
  BOOST_CHECK_EQUAL(tst_obj.ts_nfinelayers(), 0);
  BOOST_CHECK_EQUAL(tst_obj.ts_n_thermal_coeffs(), 0);
  BOOST_CHECK_CLOSE(tst_obj.ts_vlidort_accuracy(), 0, 1e-10);
  BOOST_CHECK_CLOSE(tst_obj.ts_asymtx_tolerance(), 0, 1e-10);
  BOOST_CHECK_EQUAL(tst_obj.ts_nlayers_noms(), 0);
  BOOST_CHECK_EQUAL(tst_obj.ts_nlayers_cutoff(), 0);
  BOOST_CHECK_EQUAL(tst_obj.ts_tf_maxiter(), 0);
  BOOST_CHECK_CLOSE(tst_obj.ts_tf_criterion(), 0, 1e-10);
  BOOST_CHECK_CLOSE(tst_obj.ts_toa_illumination(), 0, 1e-10);
  BOOST_CHECK_CLOSE(tst_obj.ts_boa_illumination(), 0, 1e-10);
  
}

BOOST_AUTO_TEST_CASE(vlidort_fixed_sunrays)
{
  // Used for checking dimensions
  VLidort_Pars lid_pars = VLidort_Pars::instance();

  // Test constructor
  VLidort_Fixed_Sunrays tst_obj = VLidort_Fixed_Sunrays();

  // Test memory space sizes
  tst_obj.check_byte_sizes(); // throws exception on mismatch

  // Test variable shapes
  

  // Test initialization
  BOOST_CHECK_CLOSE(tst_obj.ts_flux_factor(), 0, 1e-10);
  
}

BOOST_AUTO_TEST_CASE(vlidort_fixed_uservalues)
{
  // Used for checking dimensions
  VLidort_Pars lid_pars = VLidort_Pars::instance();

  // Test constructor
  VLidort_Fixed_Uservalues tst_obj = VLidort_Fixed_Uservalues();

  // Test memory space sizes
  tst_obj.check_byte_sizes(); // throws exception on mismatch

  // Test variable shapes
  

  // Test initialization
  BOOST_CHECK_EQUAL(tst_obj.ts_n_user_levels(), 0);
  
}

BOOST_AUTO_TEST_CASE(vlidort_fixed_chapman)
{
  // Used for checking dimensions
  VLidort_Pars lid_pars = VLidort_Pars::instance();

  // Test constructor
  VLidort_Fixed_Chapman tst_obj = VLidort_Fixed_Chapman();

  // Test memory space sizes
  tst_obj.check_byte_sizes(); // throws exception on mismatch

  // Test variable shapes
  BOOST_CHECK_EQUAL(tst_obj.ts_height_grid().extent(0), lid_pars.maxlayers+1);
  BOOST_CHECK_EQUAL(tst_obj.ts_pressure_grid().extent(0), lid_pars.maxlayers+1);
  BOOST_CHECK_EQUAL(tst_obj.ts_temperature_grid().extent(0), lid_pars.maxlayers+1);
  BOOST_CHECK_EQUAL(tst_obj.ts_finegrid().extent(0), lid_pars.maxlayers);
  

  // Test initialization
  blitz::Array<double, 1> ts_height_grid_exp(tst_obj.ts_height_grid().shape());
  ts_height_grid_exp = 0;
  BOOST_CHECK_MATRIX_CLOSE_TOL(tst_obj.ts_height_grid(), ts_height_grid_exp, 1e-10);
  blitz::Array<double, 1> ts_pressure_grid_exp(tst_obj.ts_pressure_grid().shape());
  ts_pressure_grid_exp = 0;
  BOOST_CHECK_MATRIX_CLOSE_TOL(tst_obj.ts_pressure_grid(), ts_pressure_grid_exp, 1e-10);
  blitz::Array<double, 1> ts_temperature_grid_exp(tst_obj.ts_temperature_grid().shape());
  ts_temperature_grid_exp = 0;
  BOOST_CHECK_MATRIX_CLOSE_TOL(tst_obj.ts_temperature_grid(), ts_temperature_grid_exp, 1e-10);
  blitz::Array<int, 1> ts_finegrid_exp(tst_obj.ts_finegrid().shape());
  ts_finegrid_exp = 0;
  BOOST_CHECK_MATRIX_CLOSE_TOL(tst_obj.ts_finegrid(), ts_finegrid_exp, 1e-10);
  BOOST_CHECK_CLOSE(tst_obj.ts_rfindex_parameter(), 0, 1e-10);
  
}

BOOST_AUTO_TEST_CASE(vlidort_fixed_optical)
{
  // Used for checking dimensions
  VLidort_Pars lid_pars = VLidort_Pars::instance();

  // Test constructor
  VLidort_Fixed_Optical tst_obj = VLidort_Fixed_Optical();

  // Test memory space sizes
  tst_obj.check_byte_sizes(); // throws exception on mismatch

  // Test variable shapes
  BOOST_CHECK_EQUAL(tst_obj.ts_deltau_vert_input().extent(0), lid_pars.maxlayers);
  BOOST_CHECK_EQUAL(tst_obj.ts_greekmat_total_input().extent(0), lid_pars.maxmoments_input+1);
  BOOST_CHECK_EQUAL(tst_obj.ts_greekmat_total_input().extent(1), lid_pars.maxlayers);
  BOOST_CHECK_EQUAL(tst_obj.ts_greekmat_total_input().extent(2), lid_pars.maxstokes_sq);
  BOOST_CHECK_EQUAL(tst_obj.ts_fmatrix_up().extent(0), lid_pars.maxlayers);
  BOOST_CHECK_EQUAL(tst_obj.ts_fmatrix_up().extent(1), lid_pars.max_geometries);
  BOOST_CHECK_EQUAL(tst_obj.ts_fmatrix_up().extent(2), lid_pars.6);
  BOOST_CHECK_EQUAL(tst_obj.ts_fmatrix_dn().extent(0), lid_pars.maxlayers);
  BOOST_CHECK_EQUAL(tst_obj.ts_fmatrix_dn().extent(1), lid_pars.max_geometries);
  BOOST_CHECK_EQUAL(tst_obj.ts_fmatrix_dn().extent(2), lid_pars.6);
  BOOST_CHECK_EQUAL(tst_obj.ts_thermal_bb_input().extent(0), lid_pars.maxlayers+1);
  

  // Test initialization
  blitz::Array<double, 1> ts_deltau_vert_input_exp(tst_obj.ts_deltau_vert_input().shape());
  ts_deltau_vert_input_exp = 0;
  BOOST_CHECK_MATRIX_CLOSE_TOL(tst_obj.ts_deltau_vert_input(), ts_deltau_vert_input_exp, 1e-10);
  blitz::Array<double, 3> ts_greekmat_total_input_exp(tst_obj.ts_greekmat_total_input().shape());
  ts_greekmat_total_input_exp = 0;
  BOOST_CHECK_MATRIX_CLOSE_TOL(tst_obj.ts_greekmat_total_input(), ts_greekmat_total_input_exp, 1e-10);
  blitz::Array<double, 3> ts_fmatrix_up_exp(tst_obj.ts_fmatrix_up().shape());
  ts_fmatrix_up_exp = 0;
  BOOST_CHECK_MATRIX_CLOSE_TOL(tst_obj.ts_fmatrix_up(), ts_fmatrix_up_exp, 1e-10);
  blitz::Array<double, 3> ts_fmatrix_dn_exp(tst_obj.ts_fmatrix_dn().shape());
  ts_fmatrix_dn_exp = 0;
  BOOST_CHECK_MATRIX_CLOSE_TOL(tst_obj.ts_fmatrix_dn(), ts_fmatrix_dn_exp, 1e-10);
  BOOST_CHECK_CLOSE(tst_obj.ts_lambertian_albedo(), 0, 1e-10);
  blitz::Array<double, 1> ts_thermal_bb_input_exp(tst_obj.ts_thermal_bb_input().shape());
  ts_thermal_bb_input_exp = 0;
  BOOST_CHECK_MATRIX_CLOSE_TOL(tst_obj.ts_thermal_bb_input(), ts_thermal_bb_input_exp, 1e-10);
  BOOST_CHECK_CLOSE(tst_obj.ts_surface_bb_input(), 0, 1e-10);
  BOOST_CHECK_CLOSE(tst_obj.ts_atmos_wavelength(), 0, 1e-10);
  
}

BOOST_AUTO_TEST_CASE(vlidort_fixed_write)
{
  // Used for checking dimensions
  VLidort_Pars lid_pars = VLidort_Pars::instance();

  // Test constructor
  VLidort_Fixed_Write tst_obj = VLidort_Fixed_Write();

  // Test memory space sizes
  tst_obj.check_byte_sizes(); // throws exception on mismatch

  // Test variable shapes
  

  // Test initialization
  BOOST_CHECK_EQUAL(tst_obj.ts_do_debug_write(), false);
  BOOST_CHECK_EQUAL(tst_obj.ts_do_write_input(), false);
  BOOST_CHECK_EQUAL(tst_obj.ts_do_write_scenario(), false);
  BOOST_CHECK_EQUAL(tst_obj.ts_do_write_fourier(), false);
  BOOST_CHECK_EQUAL(tst_obj.ts_do_write_results(), false);
  
}

BOOST_AUTO_TEST_CASE(vlidort_fixed_inputs)
{
  // Used for checking dimensions
  VLidort_Pars lid_pars = VLidort_Pars::instance();

  // Test constructor
  VLidort_Fixed_Inputs tst_obj = VLidort_Fixed_Inputs();

  // Test memory space sizes
  tst_obj.check_byte_sizes(); // throws exception on mismatch

  // Test variable shapes
  

  // Test initialization
  
}

BOOST_AUTO_TEST_CASE(vlidort_modified_boolean)
{
  // Used for checking dimensions
  VLidort_Pars lid_pars = VLidort_Pars::instance();

  // Test constructor
  VLidort_Modified_Boolean tst_obj = VLidort_Modified_Boolean();

  // Test memory space sizes
  tst_obj.check_byte_sizes(); // throws exception on mismatch

  // Test variable shapes
  

  // Test initialization
  BOOST_CHECK_EQUAL(tst_obj.ts_do_focorr(), false);
  BOOST_CHECK_EQUAL(tst_obj.ts_do_focorr_external(), false);
  BOOST_CHECK_EQUAL(tst_obj.ts_do_focorr_nadir(), false);
  BOOST_CHECK_EQUAL(tst_obj.ts_do_focorr_outgoing(), false);
  BOOST_CHECK_EQUAL(tst_obj.ts_do_sscorr_truncation(), false);
  BOOST_CHECK_EQUAL(tst_obj.ts_do_sscorr_usefmat(), false);
  BOOST_CHECK_EQUAL(tst_obj.ts_do_external_wleave(), false);
  BOOST_CHECK_EQUAL(tst_obj.ts_do_double_convtest(), false);
  BOOST_CHECK_EQUAL(tst_obj.ts_do_solar_sources(), false);
  BOOST_CHECK_EQUAL(tst_obj.ts_do_classical_solution(), false);
  BOOST_CHECK_EQUAL(tst_obj.ts_do_refractive_geometry(), false);
  BOOST_CHECK_EQUAL(tst_obj.ts_do_chapman_function(), false);
  BOOST_CHECK_EQUAL(tst_obj.ts_do_rayleigh_only(), false);
  BOOST_CHECK_EQUAL(tst_obj.ts_do_deltam_scaling(), false);
  BOOST_CHECK_EQUAL(tst_obj.ts_do_solution_saving(), false);
  BOOST_CHECK_EQUAL(tst_obj.ts_do_bvp_telescoping(), false);
  BOOST_CHECK_EQUAL(tst_obj.ts_do_user_vzangles(), false);
  BOOST_CHECK_EQUAL(tst_obj.ts_do_additional_mvout(), false);
  BOOST_CHECK_EQUAL(tst_obj.ts_do_mvout_only(), false);
  BOOST_CHECK_EQUAL(tst_obj.ts_do_thermal_transonly(), false);
  BOOST_CHECK_EQUAL(tst_obj.ts_do_observation_geometry(), false);
  BOOST_CHECK_EQUAL(tst_obj.ts_do_doublet_geometry(), false);
  
}

BOOST_AUTO_TEST_CASE(vlidort_modified_control)
{
  // Used for checking dimensions
  VLidort_Pars lid_pars = VLidort_Pars::instance();

  // Test constructor
  VLidort_Modified_Control tst_obj = VLidort_Modified_Control();

  // Test memory space sizes
  tst_obj.check_byte_sizes(); // throws exception on mismatch

  // Test variable shapes
  

  // Test initialization
  BOOST_CHECK_EQUAL(tst_obj.ts_ngreek_moments_input(), 0);
  
}

BOOST_AUTO_TEST_CASE(vlidort_modified_sunrays)
{
  // Used for checking dimensions
  VLidort_Pars lid_pars = VLidort_Pars::instance();

  // Test constructor
  VLidort_Modified_Sunrays tst_obj = VLidort_Modified_Sunrays();

  // Test memory space sizes
  tst_obj.check_byte_sizes(); // throws exception on mismatch

  // Test variable shapes
  BOOST_CHECK_EQUAL(tst_obj.ts_szangles().extent(0), lid_pars.max_szangles);
  

  // Test initialization
  BOOST_CHECK_EQUAL(tst_obj.ts_n_szangles(), 0);
  blitz::Array<double, 1> ts_szangles_exp(tst_obj.ts_szangles().shape());
  ts_szangles_exp = 0;
  BOOST_CHECK_MATRIX_CLOSE_TOL(tst_obj.ts_szangles(), ts_szangles_exp, 1e-10);
  
}

BOOST_AUTO_TEST_CASE(vlidort_modified_uservalues)
{
  // Used for checking dimensions
  VLidort_Pars lid_pars = VLidort_Pars::instance();

  // Test constructor
  VLidort_Modified_Uservalues tst_obj = VLidort_Modified_Uservalues();

  // Test memory space sizes
  tst_obj.check_byte_sizes(); // throws exception on mismatch

  // Test variable shapes
  BOOST_CHECK_EQUAL(tst_obj.ts_user_relazms().extent(0), lid_pars.max_user_relazms);
  BOOST_CHECK_EQUAL(tst_obj.ts_user_vzangles_input().extent(0), lid_pars.max_user_vzangles);
  BOOST_CHECK_EQUAL(tst_obj.ts_user_levels().extent(0), lid_pars.max_user_levels);
  BOOST_CHECK_EQUAL(tst_obj.ts_user_obsgeoms_input().extent(0), lid_pars.max_user_obsgeoms);
  BOOST_CHECK_EQUAL(tst_obj.ts_user_obsgeoms_input().extent(1), lid_pars.3);
  BOOST_CHECK_EQUAL(tst_obj.ts_user_doublets().extent(0), lid_pars.max_user_vzangles);
  BOOST_CHECK_EQUAL(tst_obj.ts_user_doublets().extent(1), lid_pars.2);
  

  // Test initialization
  BOOST_CHECK_EQUAL(tst_obj.ts_n_user_relazms(), 0);
  blitz::Array<double, 1> ts_user_relazms_exp(tst_obj.ts_user_relazms().shape());
  ts_user_relazms_exp = 0;
  BOOST_CHECK_MATRIX_CLOSE_TOL(tst_obj.ts_user_relazms(), ts_user_relazms_exp, 1e-10);
  BOOST_CHECK_EQUAL(tst_obj.ts_n_user_vzangles(), 0);
  blitz::Array<double, 1> ts_user_vzangles_input_exp(tst_obj.ts_user_vzangles_input().shape());
  ts_user_vzangles_input_exp = 0;
  BOOST_CHECK_MATRIX_CLOSE_TOL(tst_obj.ts_user_vzangles_input(), ts_user_vzangles_input_exp, 1e-10);
  blitz::Array<double, 1> ts_user_levels_exp(tst_obj.ts_user_levels().shape());
  ts_user_levels_exp = 0;
  BOOST_CHECK_MATRIX_CLOSE_TOL(tst_obj.ts_user_levels(), ts_user_levels_exp, 1e-10);
  BOOST_CHECK_CLOSE(tst_obj.ts_geometry_specheight(), 0, 1e-10);
  BOOST_CHECK_EQUAL(tst_obj.ts_n_user_obsgeoms(), 0);
  blitz::Array<double, 2> ts_user_obsgeoms_input_exp(tst_obj.ts_user_obsgeoms_input().shape());
  ts_user_obsgeoms_input_exp = 0;
  BOOST_CHECK_MATRIX_CLOSE_TOL(tst_obj.ts_user_obsgeoms_input(), ts_user_obsgeoms_input_exp, 1e-10);
  BOOST_CHECK_EQUAL(tst_obj.ts_n_user_doublets(), 0);
  blitz::Array<double, 2> ts_user_doublets_exp(tst_obj.ts_user_doublets().shape());
  ts_user_doublets_exp = 0;
  BOOST_CHECK_MATRIX_CLOSE_TOL(tst_obj.ts_user_doublets(), ts_user_doublets_exp, 1e-10);
  
}

BOOST_AUTO_TEST_CASE(vlidort_modified_chapman)
{
  // Used for checking dimensions
  VLidort_Pars lid_pars = VLidort_Pars::instance();

  // Test constructor
  VLidort_Modified_Chapman tst_obj = VLidort_Modified_Chapman();

  // Test memory space sizes
  tst_obj.check_byte_sizes(); // throws exception on mismatch

  // Test variable shapes
  

  // Test initialization
  BOOST_CHECK_CLOSE(tst_obj.ts_earth_radius(), 0, 1e-10);
  
}

BOOST_AUTO_TEST_CASE(vlidort_modified_optical)
{
  // Used for checking dimensions
  VLidort_Pars lid_pars = VLidort_Pars::instance();

  // Test constructor
  VLidort_Modified_Optical tst_obj = VLidort_Modified_Optical();

  // Test memory space sizes
  tst_obj.check_byte_sizes(); // throws exception on mismatch

  // Test variable shapes
  BOOST_CHECK_EQUAL(tst_obj.ts_omega_total_input().extent(0), lid_pars.maxlayers);
  

  // Test initialization
  blitz::Array<double, 1> ts_omega_total_input_exp(tst_obj.ts_omega_total_input().shape());
  ts_omega_total_input_exp = 0;
  BOOST_CHECK_MATRIX_CLOSE_TOL(tst_obj.ts_omega_total_input(), ts_omega_total_input_exp, 1e-10);
  
}

BOOST_AUTO_TEST_CASE(vlidort_modified_inputs)
{
  // Used for checking dimensions
  VLidort_Pars lid_pars = VLidort_Pars::instance();

  // Test constructor
  VLidort_Modified_Inputs tst_obj = VLidort_Modified_Inputs();

  // Test memory space sizes
  tst_obj.check_byte_sizes(); // throws exception on mismatch

  // Test variable shapes
  

  // Test initialization
  
}



BOOST_AUTO_TEST_SUITE_END()
       
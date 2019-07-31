#include "first_order_driver.h"
#include "unit_test_support.h"

//#include "ground.h"
#include "spurr_brdf_types.h"
#include "planck.h"
#include "lidort_driver.h"
#include "old_constant.h"

using namespace FullPhysics;
using namespace blitz;

BOOST_FIXTURE_TEST_SUITE(first_order_driver, GlobalFixture)

void test_first_order(int surface_type, ArrayAd<double, 1>& surface_params, ArrayAd<double, 1>& taug, ArrayAd<double, 1>& taur, Array<double, 1>& pert_atm, Array<double, 1>& pert_surf, bool do_solar, bool do_thermal, bool debug_output)
{
  blitz::Array<double, 1> sza, zen, azm;
  sza.resize(1); zen.resize(1); azm.resize(1);
  sza = 0.0;
  zen = 0.0;
  azm = 0.0;
  bool pure_nadir = false;
  int nstreams = 1;
  int nmoms = 2;

  int nlayer = taug.rows();
  int nparam = taug.number_variable();

  // Set up convenience ranges
  Range all = Range::all();
  Range rjac(0,nparam-1);
  Range rlay(0,nlayer-1);
  Range rsurf(0,surface_params.number_variable()-1);
  double refl_fo;
  double refl_lid;

  blitz::Array<double, 2> jac_atm_fo;
  blitz::Array<double, 1> jac_surf_param_fo;
  double jac_surf_temp_fo;
  blitz::Array<double, 1> jac_atm_temp_fo;

  blitz::Array<double, 2> jac_atm_lid;
  blitz::Array<double, 1> jac_surf_param_lid;
  double jac_surf_temp_lid;
  blitz::Array<double, 1> jac_atm_temp_lid;

  FirstOrderDriver fo_driver = FirstOrderDriver(nlayer, surface_type, nstreams, nmoms);

  // Use plane parallel since it has to be off from LIDORT
  fo_driver.set_plane_parallel();

  // Use LIDORT for comparison
  bool do_multiple_scattering_only = false;
  LidortRtDriver lidort_driver = LidortRtDriver(nstreams, nmoms, 
                                                do_multiple_scattering_only,
                                                surface_type, zen, pure_nadir,
                                                do_solar, do_thermal);  

  // Plane-parallel
  lidort_driver.set_plane_parallel();

  // Set up LIDORT to turn on single scattering calculations
  auto fbool = lidort_driver.lidort_interface()->lidort_fixin().f_bool();
  auto mbool = lidort_driver.lidort_interface()->lidort_modin().mbool();
  fbool.ts_do_ssfull(true);
  mbool.ts_do_sscorr_nadir(true);
  mbool.ts_do_deltam_scaling(false);

  // Simple height grid evenly spaced
  Array<double, 1> heights(nlayer+1);
  heights(0) = 100;
  for(int hidx = 1; hidx < nlayer+1; hidx++) {
    heights(hidx) = heights(hidx-1) - heights(0)/nlayer;
  }

  // Set up atmosphere jacobians to be one for taug and the other for taur
  taug.jacobian()(all, 0) = 1;
  taur.jacobian()(all, 1) = 1;

  ArrayAd<double, 1> od(nlayer, nparam);
  ArrayAd<double, 1> ssa(nlayer, nparam);

  for(int lay_idx = 0; lay_idx < nlayer; lay_idx++) {
    od(lay_idx) = taur(lay_idx) + taug(lay_idx);
    ssa(lay_idx) = taur(lay_idx) / od(lay_idx);
  }
  
  // No aerosols, and depolization factor = 0 
  // so simplified phase function moments:
  ArrayAd<double, 2> pf(3, nlayer, nparam);
  pf = 0.0;
  pf(0, all) = 1.0;
  pf(2, all) = 0.5;

  // Set up thermal inputs if enabled
  double bb_surface = 0.0;
  Array<double, 1> bb_atm(nlayer + 1);
  if (do_thermal) {
      // Just some nominal values to calculate the planck function
      double wn = 568.69;
      double temperature = 290.0;
      bb_surface = planck(wn, temperature);
      bb_atm = planck(wn, temperature);
  }
  
  // Set up jacobians
  surface_params.jacobian() = 1.0;

  // Make copies of surface params as it will be modified
  // internally by the routines
  ArrayAd<double, 1> lid_surface_params;
  lid_surface_params = surface_params;
  ArrayAd<double, 1> fo_surface_params;
  fo_surface_params = surface_params;

  // Run lidort and first order to generate values for comparison
  lidort_driver.reflectance_and_jacobian_calculate(heights, sza(0), zen(0), azm(0),
                                                   surface_type, lid_surface_params,
                                                   od, ssa, pf, refl_lid,
                                                   jac_atm_lid, jac_surf_param_lid, jac_surf_temp_lid, jac_atm_temp_lid,
                                                   bb_surface, bb_atm);
  fo_driver.reflectance_and_jacobian_calculate(heights, sza(0), zen(0), azm(0),
                                               surface_type, fo_surface_params,
                                               od, ssa, pf, refl_fo,
                                               jac_atm_fo, jac_surf_param_fo, jac_surf_temp_fo, jac_atm_temp_fo,
                                               bb_surface, bb_atm);

  if(debug_output) {
    std::cerr << "refl_lid = " << refl_lid << std::endl
              << "refl_fo  = " << refl_fo << std::endl;

  }

  BOOST_CHECK_CLOSE(refl_lid, refl_fo, 5e-3);

  blitz::Array<double, 2> jac_atm_fd(nparam, nlayer);
  double refl_fd;

  for(int l_idx = 0; l_idx < nlayer; l_idx++) {
    for(int p_idx = 0; p_idx < nparam; p_idx++) {
      blitz::Array<double,1> taug_pert( taug.value().copy() );
      blitz::Array<double,1> taur_pert( taur.value().copy() );

      blitz::Array<double,1> od_pert( od.value().shape() );
      blitz::Array<double,1> ssa_pert( ssa.value().shape() );
      blitz::Array<double,2> pf_pert( pf.value().copy() );
   
      switch (p_idx) {
      case 0:
        taug_pert(l_idx) += pert_atm(p_idx);
        break;
      case 1:
        taur_pert(l_idx) += pert_atm(p_idx);
        break;
      }

      od_pert = taur_pert + taug_pert;
      ssa_pert = taur_pert / od_pert;

      refl_fd = fo_driver.reflectance_calculate(heights, sza(0), zen(0), azm(0),
                                                surface_type, surface_params.value().copy(),
                                                od_pert, ssa_pert, pf_pert,
                                                bb_surface, bb_atm);

      jac_atm_fd(p_idx, l_idx) = (refl_fd - refl_fo) / pert_atm(p_idx);
    }
  }

  if(debug_output) {
    std::cerr << setprecision(8)
              << "jac_atm_lid = " << jac_atm_lid(rjac,rlay).transpose(1,0)
              << "jac_atm_fd = " << jac_atm_fd(rjac, rlay).transpose(1,0)
              << "jac_atm_fo = " << jac_atm_fo(rjac,rlay).transpose(1,0) << std::endl;
  }

}

void test_first_order_lambertian(bool do_solar, bool do_thermal, bool debug_output)
{
  int nlayer = 2;
  int nparam = 2;
  ArrayAd<double, 1> surface_params(1, 1);
  ArrayAd<double, 1> taug(nlayer, nparam);
  ArrayAd<double, 1> taur(nlayer, nparam);

  // Use simple lambertian throughout
  int surface_type = LAMBERTIAN;

  // Check jacobian using finite derivatives
  // od, ssa
  blitz::Array<double, 1> pert_atm(2);

  // Perturbation value for surface fd checking
  Array<double, 1> pert_surf(1);

  ////////////////
  // Surface only
  if (do_thermal && !do_solar) {
      // Make emissivity a large number
      surface_params.value() = 1.0e-6;
  } else {
      surface_params.value() = 1.0;
  }

  taur = 1.0e-6/nlayer;
  taug = 1.0e-6/nlayer;

  if (debug_output) std::cerr << "Surface only" << std::endl << "----------------------------" << std::endl;  
  pert_atm = 1e-4, 1e-4;
  pert_surf = 1e-8;
  test_first_order(surface_type, surface_params, taug, taur, pert_atm, pert_surf, do_solar, do_thermal, debug_output);

  ////////////////
  // Rayleigh only

  if (do_thermal && !do_solar) {
      // For thermal mode make the surface makes emissivity close to 0. Finite difference goes screwy if this is set == 1.0
      surface_params = 0.98;
  } else {
      // 2stream will divide by 0 in thermal mode if surface param is set to identically 0 for lambertian mode
      surface_params = 1.0e-6;
  }

  taur = 2.0e-2/nlayer;
  taug = 1.0e-6/nlayer;

  if (debug_output) std::cerr << "Rayleigh only" << std::endl << "----------------------------" << std::endl;  
  pert_atm = 1e-3, -1e-4;
  pert_surf = 1e-8;
  test_first_order(surface_type, surface_params, taug, taur, pert_atm, pert_surf, do_solar, do_thermal, debug_output);

  ////////////////
  // Gas + Surface
  if (do_thermal && !do_solar) {
      surface_params = 1.0e-6;
  } else {
      surface_params = 1.0;
  }

  taur = 1.0e-6/nlayer;
  taug = 1.0/nlayer;

  if (debug_output) std::cerr << "Gas + Surface" << std::endl << "----------------------------" << std::endl;  
  pert_atm = 1e-4, 1e-4;
  pert_surf = 1e-8;
  test_first_order(surface_type, surface_params, taug, taur, pert_atm, pert_surf, do_solar, do_thermal, debug_output);

}

BOOST_AUTO_TEST_CASE(lambertian_solar)
{
  bool do_solar = true;
  bool do_thermal = false;
  bool debug_output = true;

  test_first_order_lambertian(do_solar, do_thermal, debug_output);
}

BOOST_AUTO_TEST_SUITE_END()

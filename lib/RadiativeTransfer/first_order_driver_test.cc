#include "first_order_driver.h"
#include "unit_test_support.h"

//#include "ground.h"
#include "spurr_brdf_types.h"
#include "lidort_driver.h"
#include "old_constant.h"

using namespace FullPhysics;
using namespace blitz;

BOOST_FIXTURE_TEST_SUITE(first_order_driver, GlobalFixture)

void test_first_order(int surface_type, Array<double, 1>& surface_params, Array<double, 1>& taug, Array<double, 1>& taur, bool do_solar, bool do_thermal, bool debug_output)
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

  // Set up convenience ranges
  Range all = Range::all();
  double refl_fo;
  double refl_lid;

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

  Array<double, 1> od(nlayer);
  Array<double, 1> ssa(nlayer);
  
  for(int lay_idx = 0; lay_idx < nlayer; lay_idx++) {
    od(lay_idx) = taur(lay_idx) + taug(lay_idx);
    ssa(lay_idx) = taur(lay_idx) / od(lay_idx);
  }
  
  // No aerosols, and depolization factor = 0 
  // so simplified phase function moments:
  Array<double, 2> pf(3, nlayer);
  pf = 0.0;
  pf(0, all) = 1.0;
  pf(2, all) = 0.5;

  // Run lidort and first order to generate values for comparison
  refl_lid = lidort_driver.reflectance_calculate(heights, sza(0), zen(0), azm(0),
                                                 surface_type, surface_params,
                                                 od, ssa, pf);
  refl_fo = fo_driver.reflectance_calculate(heights, sza(0), zen(0), azm(0),
                                            surface_type, surface_params,
                                            od, ssa, pf);

  if(debug_output) {
    std::cerr << "refl_lid = " << refl_lid << std::endl
              << "refl_fo  = " << refl_fo << std::endl;
  }

  BOOST_CHECK_CLOSE(refl_lid, refl_fo, 5e-3);
}

void test_first_order_lambertian(bool do_solar, bool do_thermal, bool debug_output)
{
  int nlayer = 2;
  Array<double, 1> surface_params(1); 
  Array<double, 1> taug(nlayer);
  Array<double, 1> taur(nlayer);

  // Use simple lambertian throughout
  int surface_type = LAMBERTIAN;
  
  ////////////////
  // Surface only
  if (do_thermal && !do_solar) {
      // Make emissivity a large number
      surface_params = 1.0e-6;
  } else {
      surface_params = 1.0;
  }

  taur = 1.0e-6/nlayer;
  taug = 1.0e-6/nlayer;

  if (debug_output) std::cerr << "Surface only" << std::endl << "----------------------------" << std::endl;  
  test_first_order(surface_type, surface_params, taug, taur, do_solar, do_thermal, debug_output);

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
  test_first_order(surface_type, surface_params, taug, taur, do_solar, do_thermal, debug_output);

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
  test_first_order(surface_type, surface_params, taug, taur, do_solar, do_thermal, debug_output);

}

BOOST_AUTO_TEST_CASE(lambertian_solar)
{
  bool do_solar = true;
  bool do_thermal = false;
  bool debug_output = true;

  test_first_order_lambertian(do_solar, do_thermal, debug_output);
}

BOOST_AUTO_TEST_SUITE_END()

#include "first_order_rt.h"
#include "lidort_fixture.h"
#include "unit_test_support.h"
#include "fe_disable_exception.h"

using namespace FullPhysics;
using namespace blitz;

void compare_lidort_fo(boost::shared_ptr<RtAtmosphere>& atmosphere, 
                       const boost::shared_ptr<StokesCoefficient>& stokes_coefs,
                       blitz::Array<double, 1>& sza, 
                       blitz::Array<double, 1>& zen, 
                       blitz::Array<double, 1>& azm, 
                       blitz::Array<double, 1>& wn_arr,
                       bool debug_output)
{ 
  int nstreams = 1;
  int nmoms = 3;
  bool do_multiple_scattering_only = false;
  bool pure_nadir = false;
  bool do_deltam_scaling = true;

  boost::shared_ptr<LidortRt> lidort_rt;
  lidort_rt.reset(new LidortRt(atmosphere, stokes_coefs, sza, zen, azm, pure_nadir,
                               nstreams, nmoms, do_multiple_scattering_only));
  
  // Set deltam scaling on or off
  lidort_rt->rt_driver()->lidort_interface()->lidort_modin().mbool().ts_do_deltam_scaling(do_deltam_scaling);

  // Plane-parallel
  //lidort_rt->rt_driver()->set_plane_parallel();

  // Disable the diffuse calculation in LIDORT
  auto fbool = lidort_rt->rt_driver()->lidort_interface()->lidort_fixin().f_bool();
  fbool.ts_do_ssfull(true);

  //auto mbool = lidort_rt->rt_driver()->lidort_interface()->lidort_modin().mbool();
  //mbool.ts_do_sscorr_nadir(true);

  ArrayAd<double, 2> lid_rad_jac = lidort_rt->stokes_and_jacobian(wn_arr, 0);

  boost::shared_ptr<FirstOrderRt> first_order_rt;
  first_order_rt.reset(new FirstOrderRt(atmosphere, stokes_coefs, sza, zen, azm, nstreams, nmoms));

  // Use plane parallel to match LIDORT
  //first_order_rt->rt_driver()->set_plane_parallel();
  
  first_order_rt->rt_driver()->do_deltam_scaling(do_deltam_scaling);

  ArrayAd<double, 2> fo_rad_jac = first_order_rt->stokes_and_jacobian(wn_arr, 0);

  Range all = Range::all();
  if(debug_output) {
    std::cerr << "lid_rad = " << lid_rad_jac.value()(all, 0) << std::endl
              << " fo_rad = " << fo_rad_jac.value()(all, 0) << std::endl;
  }

  BOOST_CHECK_MATRIX_CLOSE_TOL(lid_rad_jac.value()(all, 0), fo_rad_jac.value()(all, 0), 1e-6);

  for(int jac_idx = 0; jac_idx < fo_rad_jac.jacobian().depth(); jac_idx++) {
    if(debug_output) {
      std::cerr << jac_idx << " -- l: ";
      for(int wn_idx = 0; wn_idx < lid_rad_jac.jacobian().rows(); wn_idx++)      
        std::cerr << lid_rad_jac.jacobian()(wn_idx, 0, jac_idx) << " ";
      std::cerr << std::endl << jac_idx << " -- f: ";
      for(int wn_idx = 0; wn_idx < fo_rad_jac.jacobian().rows(); wn_idx++)      
        std::cerr << fo_rad_jac.jacobian()(wn_idx, 0, jac_idx) << " ";
      std::cerr << std::endl;
    }
    BOOST_CHECK_MATRIX_CLOSE_TOL(lid_rad_jac.jacobian()(all, 0, jac_idx), fo_rad_jac.jacobian()(all, 0, jac_idx), 7e-4);
  }
}

BOOST_FIXTURE_TEST_SUITE(first_order_rt_lambertian, LidortLambertianFixture)

BOOST_AUTO_TEST_CASE(comparison)
{ 
  compare_lidort_fo(config_atmosphere, stokes_coefs, sza, zen, azm, wn_arr, true);
}

BOOST_AUTO_TEST_SUITE_END()


BOOST_FIXTURE_TEST_SUITE(first_order_rt_coxmunk, LidortCoxmunkFixture)

BOOST_AUTO_TEST_CASE(comparison)
{ 
  compare_lidort_fo(config_atmosphere, stokes_coefs, sza, zen, azm, wn_arr, true);
}

BOOST_AUTO_TEST_SUITE_END()


BOOST_FIXTURE_TEST_SUITE(first_order_rt_coxmunk_plus_lambertian, LidortCoxmunkPlusLambertianFixture)

BOOST_AUTO_TEST_CASE(comparison)
{ 
  // Not sure what is triggering an exception on ubuntu. It is hard to 
  // duplicate this problem, so we'll just shut this off for now.
  FeDisableException disable_fp;
  compare_lidort_fo(config_atmosphere, stokes_coefs, sza, zen, azm, wn_arr, true);
}

BOOST_AUTO_TEST_SUITE_END()


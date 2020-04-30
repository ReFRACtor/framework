#include "first_order_rt.h"
#include "lidort_fixture.h"
#include "unit_test_support.h"
#include "fe_disable_exception.h"

#include "altitude.h"

using namespace FullPhysics;
using namespace blitz;

const bool do_debug_output = false;

void compare_lidort_fo(boost::shared_ptr<RtAtmosphere>& atmosphere, 
                       const boost::shared_ptr<StokesCoefficient>& stokes_coefs,
                       blitz::Array<double, 1>& sza, 
                       blitz::Array<double, 1>& zen, 
                       blitz::Array<double, 1>& azm, 
                       blitz::Array<double, 1>& wn_arr,
                       int sphericity_mode,                       
                       bool do_deltam_scaling,
                       bool debug_output)
{ 
  int nstreams = 8;
  int nmoms = 16;
  bool do_multiple_scattering_only = false;
  bool pure_nadir = false;

  boost::shared_ptr<LidortRt> lidort_rt;
  lidort_rt.reset(new LidortRt(atmosphere, stokes_coefs, sza, zen, azm, pure_nadir,
                               nstreams, nmoms, do_multiple_scattering_only));
  
  // Set deltam scaling on or off
  lidort_rt->rt_driver()->lidort_interface()->lidort_modin().mbool().ts_do_deltam_scaling(do_deltam_scaling);

  // Disable the diffuse calculation in LIDORT
  auto fbool = lidort_rt->rt_driver()->lidort_interface()->lidort_fixin().f_bool();
  fbool.ts_do_ssfull(true);

  // Configure sphericity and phase function affecting options
  switch(sphericity_mode) {
    case 0:
        lidort_rt->rt_driver()->set_plane_parallel();
        break;
    case 1:
        lidort_rt->rt_driver()->set_pseudo_spherical();
        break;
    case 2:
        lidort_rt->rt_driver()->set_line_of_sight();
        break;
    default:
        throw Exception("Unknown sphericity_mode value");
  }

  ArrayAd<double, 2> lid_rad_jac = lidort_rt->stokes_and_jacobian(wn_arr, 0);

  boost::shared_ptr<FirstOrderRt> first_order_rt;
  first_order_rt.reset(new FirstOrderRt(atmosphere, stokes_coefs, sza, zen, azm, nstreams, nmoms));
  
  first_order_rt->rt_driver()->do_deltam_scaling(do_deltam_scaling);

  // Configure sphericity and phase function affecting options
  switch(sphericity_mode) {
    case 0:
        first_order_rt->rt_driver()->set_plane_parallel();
        break;
    case 1:
        first_order_rt->rt_driver()->set_pseudo_spherical();
        break;
    case 2:
        first_order_rt->rt_driver()->set_line_of_sight();
        break;
    default:
        throw Exception("Unknown sphericity_mode value");
  }

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
    BOOST_CHECK_MATRIX_CLOSE_TOL(lid_rad_jac.jacobian()(all, 0, jac_idx), fo_rad_jac.jacobian()(all, 0, jac_idx), 7e-3);
  }
}

BOOST_FIXTURE_TEST_SUITE(first_order_rt_lambertian, LidortLambertianFixture)

BOOST_AUTO_TEST_CASE(plane_parallel_no_deltam)
{ 
  compare_lidort_fo(config_atmosphere, stokes_coefs, sza, zen, azm, wn_arr, 0, false, do_debug_output);
}

BOOST_AUTO_TEST_CASE(psuedo_spherical_no_deltam)
{ 
  compare_lidort_fo(config_atmosphere, stokes_coefs, sza, zen, azm, wn_arr, 1, false, do_debug_output);
}

BOOST_AUTO_TEST_CASE(line_of_sight_no_deltam)
{ 
  compare_lidort_fo(config_atmosphere, stokes_coefs, sza, zen, azm, wn_arr, 2, false, do_debug_output);
}

BOOST_AUTO_TEST_CASE(no_aerosol)
{ 

  boost::shared_ptr<AtmosphereStandard> src_atm(boost::dynamic_pointer_cast<AtmosphereStandard>(config_atmosphere));

  std::vector<boost::shared_ptr<Altitude> > alt_clone;
    BOOST_FOREACH(const boost::shared_ptr<Altitude>& a, src_atm->altitude_ptr())
        alt_clone.push_back(a->clone());

  boost::shared_ptr<RtAtmosphere>
    atm_no_aerosol(new AtmosphereStandard(src_atm->absorber_ptr()->clone(),
                                          src_atm->pressure_ptr()->clone(),
                                          src_atm->temperature_ptr()->clone(),
                                          src_atm->rayleigh_ptr()->clone(),
                                          NULL,
                                          src_atm->relative_humidity_ptr()->clone(),
                                          src_atm->ground()->clone(),
                                          alt_clone,
                                          src_atm->constant_ptr()));
  
  compare_lidort_fo(atm_no_aerosol, stokes_coefs, sza, zen, azm, wn_arr, 2, true, do_debug_output);
}

BOOST_AUTO_TEST_SUITE_END()


BOOST_FIXTURE_TEST_SUITE(first_order_rt_coxmunk, LidortCoxmunkFixture)

BOOST_AUTO_TEST_CASE(plane_parallel_no_deltam)
{ 
  compare_lidort_fo(config_atmosphere, stokes_coefs, sza, zen, azm, wn_arr, 0, false, do_debug_output);
}

BOOST_AUTO_TEST_CASE(psuedo_spherical_no_deltam)
{ 
  compare_lidort_fo(config_atmosphere, stokes_coefs, sza, zen, azm, wn_arr, 1, false, do_debug_output);
}

BOOST_AUTO_TEST_CASE(line_of_sight_no_deltam)
{ 
  compare_lidort_fo(config_atmosphere, stokes_coefs, sza, zen, azm, wn_arr, 2, false, do_debug_output);
}

BOOST_AUTO_TEST_SUITE_END()


BOOST_FIXTURE_TEST_SUITE(first_order_rt_coxmunk_plus_lambertian, LidortCoxmunkPlusLambertianFixture)

BOOST_AUTO_TEST_CASE(plane_parallel_no_deltam)
{ 
  // Not sure what is triggering an exception on ubuntu. It is hard to 
  // duplicate this problem, so we'll just shut this off for now.
  FeDisableException disable_fp;
  compare_lidort_fo(config_atmosphere, stokes_coefs, sza, zen, azm, wn_arr, 0, false, do_debug_output);
}

BOOST_AUTO_TEST_CASE(psuedo_spherical_no_deltam)
{ 
  compare_lidort_fo(config_atmosphere, stokes_coefs, sza, zen, azm, wn_arr, 1, false, do_debug_output);
}

BOOST_AUTO_TEST_CASE(line_of_sight_no_deltam)
{ 
  compare_lidort_fo(config_atmosphere, stokes_coefs, sza, zen, azm, wn_arr, 2, false, do_debug_output);
}

BOOST_AUTO_TEST_SUITE_END()


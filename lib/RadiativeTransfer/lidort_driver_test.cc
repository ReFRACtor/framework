#include "lidort_driver.h"
#include "unit_test_support.h"
#include "lidort_fixture.h"
#include "old_constant.h"
#include "planck.h"

#include "spurr_brdf_types.h"

using namespace FullPhysics;
using namespace blitz;

bool check_brdf_inputs(boost::shared_ptr<LidortRtDriver>& lidort_driver) {
  Lidort_Brdf_Sup_Accessories brdf_check = Lidort_Brdf_Sup_Accessories(lidort_driver->brdf_interface()->brdf_sup_in_ptr(),
      lidort_driver->lidort_interface()->lidort_fixin_ptr(),
      lidort_driver->lidort_interface()->lidort_modin_ptr());
  brdf_check.brdf_input_check();

  Lidort_Exception_Handling& brdf_check_status = brdf_check.lidort_brdfcheck_status();
  Lidort_Pars lid_pars = Lidort_Pars::instance();

  if (brdf_check_status.ts_status_inputcheck() != lid_pars.lidort_success())
    std::cerr << brdf_check_status << std::endl;

  return brdf_check_status.ts_status_inputcheck() == lid_pars.lidort_success();
}

BOOST_FIXTURE_TEST_SUITE(lidort_driver_lambertian_solar, GlobalFixture)

class Fixture1 : public LidortDriverLambertianFixture
{
public:
void run_test(boost::shared_ptr<LidortRtDriver>& ldriver) {
  double refl_calc;
  double refl_expected;

  ////////////////
  // Surface only
  surface_params(0) = 1.0;

  taur = 1.0e-6/nlayer;
  taug = 1.0e-6/nlayer;

  od = taur + taug;
  ssa.value() = taur / od.value();

  refl_calc = ldriver->reflectance_calculate(heights, sza(0), zen(0), azm(0),
                                                   surface_type, surface_params,
                                                   od.value(), ssa.value(), pf.value());
  
  // Surface only = 1/pi
  // This checks % differerence, so tol is % diff
  refl_expected = 1/OldConstant::pi * surface_params(0);
  BOOST_CHECK_CLOSE(refl_expected, refl_calc, 1e-3);

  ////////////////
  // Gas + Surface
  surface_params(0) = 1.0;

  taur = 1.0e-6/nlayer;
  taug = 1.0/nlayer;

  od = taur + taug;
  ssa.value() = taur / od.value();

  refl_calc = ldriver->reflectance_calculate(heights, sza(0), zen(0), azm(0),
                                                   surface_type, surface_params,
                                                   od.value(), ssa.value(), pf.value());

  refl_expected = 1/OldConstant::pi * exp(-1/cos(sza(0))) * exp(-1/cos(zen(0)));
  BOOST_CHECK_CLOSE(refl_expected, refl_calc, 1e-3);

  ////////////////
  // Rayleigh only
  surface_params(0) = 0.0;

  taur = 2.0e-2/nlayer;
  taug = 1.0e-6/nlayer;

  od = taur + taug;
  ssa.value() = taur / od.value();

  refl_calc = ldriver->reflectance_calculate(heights, sza(0), zen(0), azm(0),
                                                   surface_type, surface_params,
                                                   od.value(), ssa.value(), pf.value());


  // Expected value from VLIDORT, this is
  // not going to agree to infinite precision
  // Note the tolerance here is in %
  refl_expected = 2.387246757232095E-003;
  BOOST_CHECK_CLOSE(refl_expected, refl_calc, 6e-3);
}
};

BOOST_AUTO_TEST_CASE(solar_sources)
{
  Fixture1 t;
  t.run_test(t.lidort_driver);
  t.lidort_driver->brdf_interface()->write_fortran_file("test.txt");
}

BOOST_AUTO_TEST_CASE(serialization)
{
  if(!have_serialize_supported())
    return;
  Fixture1 t;
  std::string d = serialize_write_string(t.lidort_driver);
  if(false)
    std::cerr << d;
  boost::shared_ptr<LidortRtDriver> ldriver_r =
    serialize_read_string<LidortRtDriver>(d);
  t.run_test(ldriver_r);
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_FIXTURE_TEST_SUITE(lidort_driver_lambertian_thermal, GlobalFixture)

class Fixture1 : public LidortDriverLambertianThermalFixture
{
public:
void run_test(boost::shared_ptr<LidortRtDriver>& ldriver) {
  double refl_calc;
  double refl_expt;

  /////////////////////////////////////////////
  // Thermal setup

  // Just some nominal values to calculate the planck function
  double wn = 568.69;
  double temperature = 290.0;
  double bb_surface = planck(wn, temperature);
  Array<double, 1> bb_atm(nlayer + 1);
  bb_atm = planck(wn, temperature);

  ////////////////////////
  // Thermal surface only
  surface_params(0) = 0.1;

  taur = 1.0e-6/nlayer;
  taug = 1.0e-6/nlayer;

  od = taur + taug;
  ssa.value() = taur / od.value();

  refl_calc = ldriver->reflectance_calculate(heights, sza(0), zen(0), azm(0),
                                                   surface_type, surface_params,
                                                   od.value(), ssa.value(), pf.value(),
                                                   bb_surface, bb_atm);

  // Expected values calculated offline in VLIDORT
  refl_expt = 0.12476970639470732;
  BOOST_CHECK_CLOSE(refl_expt, refl_calc, 1e-3);

  ////////////////////////
  // Thermal gas + surface
  surface_params(0) = 0.1;

  taur = 1.0e-6/nlayer;
  taug = 1.0/nlayer;

  od = taur + taug;
  ssa.value() = taur / od.value();

  refl_calc = ldriver->reflectance_calculate(heights, sza(0), zen(0), azm(0),
                                                   surface_type, surface_params,
                                                   od.value(), ssa.value(), pf.value(),
                                                   bb_surface, bb_atm);

  refl_expt = 0.13751414471789294;
  BOOST_CHECK_CLOSE(refl_expt, refl_calc, 1e-3);

  /////////////////////////////////////////////
  // Thermal rayleigh only, no surface, no gas
  surface_params(0) = 0.9;

  taur = 2.0e-2/nlayer;
  taug = 1.0e-6/nlayer;

  od = taur + taug;
  ssa.value() = taur / od.value();

  refl_calc = ldriver->reflectance_calculate(heights, sza(0), zen(0), azm(0),
                                                   surface_type, surface_params,
                                                   od.value(), ssa.value(), pf.value(),
                                                   bb_surface, bb_atm);

  refl_expt = 1.3967e-002;
  BOOST_CHECK_CLOSE(refl_expt, refl_calc, 7e-3);

}
};

BOOST_AUTO_TEST_CASE(thermal_emission)
{
  Fixture1 t;
  t.run_test(t.lidort_driver);
}

BOOST_AUTO_TEST_CASE(serialization)
{
  if(!have_serialize_supported())
    return;
  Fixture1 t;
  std::string d = serialize_write_string(t.lidort_driver);
  if(false)
    std::cerr << d;
  boost::shared_ptr<LidortRtDriver> ldriver_r =
    serialize_read_string<LidortRtDriver>(d);
  t.run_test(ldriver_r);
}

BOOST_AUTO_TEST_SUITE_END()
 
///////////////////////////////////////////////////////////////////////////////////////

BOOST_FIXTURE_TEST_SUITE(lidort_driver_coxmunk_simple, GlobalFixture)

class Fixture1 : public LidortDriverCoxmunkFixture
{
public:
void run_test(boost::shared_ptr<LidortRtDriver>& ldriver) {
  double refl_calc;
  double refl_expected;

  ////////////////
  // Surface only
  surface_params(0) = 1.0e-6;
  surface_params(1) = 1.334;
  surface_params(2) = 0; // Not used
  surface_params(3) = 0; // Shadowing flag
  
  taur = 1.0e-6/nlayer;
  taug = 1.0e-6/nlayer;

  od = taur + taug;
  ssa.value() = taur / od.value();

  refl_calc = ldriver->reflectance_calculate(heights, sza(0), zen(0), azm(0),
                                                   surface_type, surface_params,
                                                   od.value(), ssa.value(), pf.value());

  BOOST_CHECK_EQUAL(check_brdf_inputs(lidort_driver), true);
  
  // Value for VLIDORT
  refl_expected = 0.54319850628416033;
  BOOST_CHECK_CLOSE(refl_expected, refl_calc, 1e-3);

}
};

BOOST_AUTO_TEST_CASE(simple)
{
  Fixture1 t;
  t.run_test(t.lidort_driver);
}

BOOST_AUTO_TEST_CASE(serialization)
{
  if(!have_serialize_supported())
    return;
  Fixture1 t;
  std::string d = serialize_write_string(t.lidort_driver);
  if(false)
    std::cerr << d;
  boost::shared_ptr<LidortRtDriver> ldriver_r =
    serialize_read_string<LidortRtDriver>(d);
  t.run_test(ldriver_r);
}

BOOST_AUTO_TEST_SUITE_END()

///////////////////////////////////////////////////////////////////////////////////////

BOOST_FIXTURE_TEST_SUITE(lidort_driver_coxmunk_plus_lamb_simple, GlobalFixture)

class Fixture1 : public LidortDriverCoxmunkFixture
{
public:
void run_test(boost::shared_ptr<LidortRtDriver>& ldriver) {
  double refl_calc;
  double refl_expt;

  ////////////////
  // Surface only
  surface_params(0) = 1.0e-6;
  surface_params(1) = 1.334;
  surface_params(2) = 0.5;
  surface_params(3) = 0; // Shadowing flag

  taur = 1.0e-6/nlayer;
  taug = 1.0e-6/nlayer;

  od = taur + taug;
  ssa = 0;
  ssa.value() = taur / od.value();
  
  //////
  // Plane-parallel
  ldriver->set_plane_parallel();

  blitz::Array<double, 2> jac_atm;
  blitz::Array<double, 1> jac_surf_param;
  double jac_surf_temp;
  blitz::Array<double, 1> jac_atm_temp;

  ArrayAd<double, 1> lidort_surface(surface_params.shape(), 1);

  lidort_surface.value() = surface_params;
  lidort_surface.jacobian() = 1.0;

  ldriver->reflectance_and_jacobian_calculate(heights, sza(0), zen(0), azm(0),
                                                    surface_type, lidort_surface,
                                                    od, ssa, pf, refl_calc, 
                                                    jac_atm, jac_surf_param, jac_surf_temp, jac_atm_temp);

  BOOST_CHECK_EQUAL(check_brdf_inputs(lidort_driver), true);

  // Adjust analytic jacobians have same meaning as fd jacobians
  jac_surf_param(0) *= lidort_surface.jacobian()(0,0);
  jac_surf_param(1) *= lidort_surface.jacobian()(1,0);

  // Value for VLIDORT
  refl_expt = 0.70235315460259928;
  BOOST_CHECK_CLOSE(refl_expt, refl_calc, 1e-3);

  // Check surface jacobians against FD

  blitz::Array<double, 1> pert_values(surface_params.extent(firstDim)-1);
  pert_values = 1e-8, 1e-8, 1e-6;

  blitz::Array<double, 1> jac_surf_param_fd( jac_surf_param.extent() );
  double refl_fd;

  jac_surf_param_fd = 0.0;
  refl_fd = 0.0;

  // First check PP mode against value just computed
  for(int p_idx = 0; p_idx < pert_values.extent(firstDim); p_idx++) {
    blitz::Array<double,1> surface_params_pert( surface_params.extent() );
    surface_params_pert = surface_params;
    surface_params_pert(p_idx) += pert_values(p_idx);

    refl_fd = ldriver->reflectance_calculate(heights, sza(0), zen(0), azm(0),
                                                   surface_type, surface_params_pert,
                                                   od.value(), ssa.value(), pf.value());

    jac_surf_param_fd(p_idx) = (refl_fd - refl_calc) / pert_values(p_idx);
  }

  BOOST_CHECK_MATRIX_CLOSE_TOL(jac_surf_param, jac_surf_param_fd, 1e-7);

  // Pseudo-spherical mode FD test
  ldriver->set_pseudo_spherical();
 
  lidort_surface.value() = surface_params;
  lidort_surface.jacobian() = 1.0;

  ldriver->reflectance_and_jacobian_calculate(heights, sza(0), zen(0), azm(0),
                                                    surface_type, lidort_surface,
                                                    od, ssa, pf, refl_calc,
                                                    jac_atm, jac_surf_param, jac_surf_temp, jac_atm_temp);

  BOOST_CHECK_EQUAL(check_brdf_inputs(lidort_driver), true);

  // Adjust analytic jacobians have same meaning as fd jacobians
  jac_surf_param(0) *= lidort_surface.jacobian()(0,0);
  jac_surf_param(1) *= lidort_surface.jacobian()(1,0);

  for(int p_idx = 0; p_idx < pert_values.extent(firstDim); p_idx++) {
    blitz::Array<double,1> surface_params_pert( surface_params.extent() );
    surface_params_pert = surface_params;
    surface_params_pert(p_idx) += pert_values(p_idx);

    refl_fd = ldriver->reflectance_calculate(heights, sza(0), zen(0), azm(0),
                                                   surface_type, surface_params_pert,
                                                   od.value(), ssa.value(), pf.value());

    jac_surf_param_fd(p_idx) = (refl_fd - refl_calc) / pert_values(p_idx);
  }

  BOOST_CHECK_MATRIX_CLOSE_TOL(jac_surf_param, jac_surf_param_fd, 1e-7);

  // Line-of-site mode FD test
  sza = sza + 1e-3; // or else risk divide by zero
  azm = azm + 1e-3; // or else risk divide by zero
  zen = zen + 1e-3; // or else risk divide by zero

  ldriver->set_line_of_sight();

  lidort_surface.value() = surface_params;
  lidort_surface.jacobian() = 1.0;
  ldriver->reflectance_and_jacobian_calculate(heights, sza(0), zen(0), azm(0),
                                                 surface_type, lidort_surface,
                                                 od, ssa, pf, refl_calc,
                                                 jac_atm, jac_surf_param, jac_surf_temp, jac_atm_temp);

  BOOST_CHECK_EQUAL(check_brdf_inputs(lidort_driver), true);

  // Adjust analytic jacobians have same meaning as fd jacobians
  jac_surf_param(0) *= lidort_surface.jacobian()(0,0);
  jac_surf_param(1) *= lidort_surface.jacobian()(1,0);

  for(int p_idx = 0; p_idx < pert_values.extent(firstDim); p_idx++) {
    blitz::Array<double,1> surface_params_pert( surface_params.extent() );
    surface_params_pert = surface_params;
    surface_params_pert(p_idx) += pert_values(p_idx);

    refl_fd = ldriver->reflectance_calculate(heights, sza(0), zen(0), azm(0),
                                                   surface_type, surface_params_pert,
                                                   od.value(), ssa.value(), pf.value());

    jac_surf_param_fd(p_idx) = (refl_fd - refl_calc) / pert_values(p_idx);
  }

  BOOST_CHECK_MATRIX_CLOSE_TOL(jac_surf_param, jac_surf_param_fd, 1e-7);

}
};

BOOST_AUTO_TEST_CASE(simple)
{
  Fixture1 t;
  t.run_test(t.lidort_driver);
}

BOOST_AUTO_TEST_CASE(serialization)
{
  if(!have_serialize_supported())
    return;
  Fixture1 t;
  std::string d = serialize_write_string(t.lidort_driver);
  if(false)
    std::cerr << d;
  boost::shared_ptr<LidortRtDriver> ldriver_r =
    serialize_read_string<LidortRtDriver>(d);
  t.run_test(ldriver_r);
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_FIXTURE_TEST_SUITE(lidort_driver_brdf_veg, GlobalFixture)

class Fixture1 : public LidortDriverBrdfVegFixture
{
public:
void run_test(boost::shared_ptr<LidortRtDriver>& ldriver) {
  double refl_calc;
  blitz::Array<double, 2> jac_atm;
  blitz::Array<double, 1> jac_surf_param;
  double jac_surf_temp;
  blitz::Array<double, 1> jac_atm_temp;

  ////////////////
  // Surface only

  surface_params(0) = 1.0; // rahman kernel factor
  surface_params(1) = 0.1; // hotspot parameter
  surface_params(2) = 0.3; // asymmetry
  surface_params(3) = 1.5; // anisotropy_parameter
  surface_params(4) = 1.0; // breon kernel factor

  ArrayAd<double, 1> lidort_surface(surface_params.shape(), 1);

  lidort_surface.value() = surface_params;
  lidort_surface.jacobian() = 1.0;

  taur = 1.0e-6/nlayer;
  taug = 1.0e-6/nlayer;

  od = taur + taug;
  ssa = 0;
  ssa.value() = taur / od.value();

  // Set sza to compare with that used in the l_rad test we compare against
  sza = 0.1;

  ldriver->reflectance_and_jacobian_calculate(heights, sza(0), zen(0), azm(0),
                                                    surface_type, lidort_surface,
                                                    od, ssa, pf, refl_calc, 
                                                    jac_atm, jac_surf_param, jac_surf_temp, jac_atm_temp);

  // Compare against an offline calculated value, or could compare against value from l_rad
  double refl_expected = 0.035435854422713485;
  BOOST_CHECK_CLOSE(refl_expected, refl_calc, 1e-3);

  // Check surface jacobians against FD

  blitz::Array<double, 1> pert_values(surface_params.rows());
  pert_values = 1e-8, 1e-8, 1e-6, 1e-6, 1e-6;

  blitz::Array<double, 1> jac_surf_param_fd( jac_surf_param.extent() );
  double refl_fd;

  jac_surf_param_fd = 0.0;
  refl_fd = 0.0;

  // First check PP mode against value just computed
  for(int p_idx = 0; p_idx < pert_values.extent(firstDim); p_idx++) {
    blitz::Array<double,1> surface_params_pert( surface_params.rows() );
    surface_params_pert = surface_params;
    surface_params_pert(p_idx) += pert_values(p_idx);

    refl_fd = ldriver->reflectance_calculate(heights, sza(0), zen(0), azm(0),
                                                   surface_type, surface_params_pert,
                                                   od.value(), ssa.value(), pf.value());

    jac_surf_param_fd(p_idx) = (refl_fd - refl_calc) / pert_values(p_idx);
  }

  BOOST_CHECK_MATRIX_CLOSE_TOL(jac_surf_param, jac_surf_param_fd, 2e-7);

}
};

BOOST_AUTO_TEST_CASE(simple)
{
  Fixture1 t;
  t.run_test(t.lidort_driver);
}

BOOST_AUTO_TEST_CASE(serialization)
{
  if(!have_serialize_supported())
    return;
  Fixture1 t;
  std::string d = serialize_write_string(t.lidort_driver);
  if(false)
    std::cerr << d;
  boost::shared_ptr<LidortRtDriver> ldriver_r =
    serialize_read_string<LidortRtDriver>(d);
  t.run_test(ldriver_r);
}

BOOST_AUTO_TEST_SUITE_END()

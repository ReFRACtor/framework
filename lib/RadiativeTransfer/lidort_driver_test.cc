#include "lidort_rt_driver.h"
#include "unit_test_support.h"
#include "lidort_fixture.h"
#include "old_constant.h"
#include "planck.h"
#include "fe_disable_exception.h"
#include "thermal_jacobian_check.h"

#include "spurr_brdf_types.h"

using namespace FullPhysics;
using namespace blitz;

bool check_brdf_inputs(boost::shared_ptr<LidortRtDriver>& lidort_driver) {
  Lidort_Brdf_Sup_Accessories brdf_check = Lidort_Brdf_Sup_Accessories(lidort_driver->lidort_brdf_interface()->brdf_sup_in_ptr(),
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
  Array<double, 1> refl_calc(1);
  Array<double, 1> refl_expected(1);

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
  BOOST_CHECK_MATRIX_CLOSE_TOL(refl_expected, refl_calc, 1e-3);

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
  BOOST_CHECK_MATRIX_CLOSE_TOL(refl_expected, refl_calc, 1e-3);

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
  BOOST_CHECK_MATRIX_CLOSE_TOL(refl_expected, refl_calc, 6e-3);
}
};

BOOST_AUTO_TEST_CASE(solar_sources)
{
  Fixture1 t;
  t.run_test(t.lidort_driver);
  t.lidort_driver->lidort_brdf_interface()->write_fortran_file("test.txt");
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
  Array<double, 1> refl_calc(1);
  Array<double, 1> refl_expt(1);

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
  BOOST_CHECK_MATRIX_CLOSE_TOL(refl_expt, refl_calc, 1e-5);

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
  BOOST_CHECK_MATRIX_CLOSE_TOL(refl_expt, refl_calc, 1e-5);

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

  refl_expt = 1.3966249402729437e-2;
  BOOST_CHECK_MATRIX_CLOSE_TOL(refl_expt, refl_calc, 1e-5);

}
};

BOOST_AUTO_TEST_CASE(thermal_emission)
{
  Fixture1 t;
  t.run_test(t.lidort_driver);
}

// Blackbody weighting-function Jacobians (jac_surf_temp = d(refl)/d(surface
// blackbody), jac_atm_temp = d(refl)/d(per-level atmosphere blackbody))
// versus finite differences. Thermal emission is linear in the blackbody
// inputs, so the finite differences are effectively exact.
BOOST_AUTO_TEST_CASE(blackbody_jacobian_fd)
{
  FeDisableException disable_fp;
  Fixture1 t;

  double wn = 568.69;
  double temperature = 290.0;
  double bb_surface = planck(wn, temperature);
  Array<double, 1> bb_atm(t.nlayer + 1);
  bb_atm = planck(wn, temperature);

  // Gas + (near-black, fully emitting) surface. The surface albedo is set
  // ~0 so its emissivity is ~1 and there is no reflection of downwelling
  // atmospheric emission back up: this isolates the surface and atmospheric
  // blackbody weighting functions so each analytic Jacobian can be compared
  // against a clean finite difference.
  double taur = 1.0e-6/t.nlayer;
  double taug = 1.0/t.nlayer;

  ArrayAd<double, 1> od(t.nlayer, 1), ssa(t.nlayer, 1);
  od.value() = taur + taug;
  od.jacobian() = 0.0;
  ssa.value() = taur / od.value();
  ssa.jacobian() = 0.0;

  ArrayAd<double, 3> pf(t.pf);
  pf.jacobian() = 0.0;

  ArrayAd<double, 1> surface(1, 1);
  surface.value()(0) = 1.0e-6;
  surface.jacobian()(0, 0) = 0.0;

  Array<double, 1> refl;
  Array<double, 3> jac_atm;
  Array<double, 2> jac_surf_param, jac_atm_temp;
  Array<double, 1> jac_surf_temp;

  t.lidort_driver->reflectance_and_jacobian_calculate(
      t.heights, t.sza(0), t.azm(0), t.zen(0), t.surface_type, surface,
      od, ssa, pf, refl, jac_atm, jac_surf_param, jac_surf_temp, jac_atm_temp,
      bb_surface, bb_atm);

  double pert = 1e-3;

  // Surface blackbody jacobian
  {
    Array<double, 1> refl_pert(t.lidort_driver->reflectance_calculate(
        t.heights, t.sza(0), t.azm(0), t.zen(0), t.surface_type, surface.value(),
        od.value(), ssa.value(), pf.value(), bb_surface + pert, bb_atm));
    double jac_fd = (refl_pert(0) - refl(0)) / pert;
    BOOST_CHECK_CLOSE(jac_surf_temp(0), jac_fd, 1e-2);
  }

  // Per-level atmosphere blackbody jacobian
  for(int lev = 0; lev < bb_atm.rows(); ++lev) {
    Array<double, 1> bb_atm_pert(bb_atm.copy());
    bb_atm_pert(lev) += pert;
    Array<double, 1> refl_pert(t.lidort_driver->reflectance_calculate(
        t.heights, t.sza(0), t.azm(0), t.zen(0), t.surface_type, surface.value(),
        od.value(), ssa.value(), pf.value(), bb_surface, bb_atm_pert));
    double jac_fd = (refl_pert(0) - refl(0)) / pert;
    BOOST_CHECK_CLOSE(jac_atm_temp(lev, 0), jac_fd, 1e-2);
  }
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

//-----------------------------------------------------------------------
// Downwelling thermal radiance at the bottom of the atmosphere (an
// up-looking observer at the surface, nadir-zenith view). Expected values
// are from the independent successive-orders solution
// (reference/rt/vlidort_thermal_test/sos_thermal_check.py). The gas case is
// also analytic: B (1 - exp(-tau)). The Rayleigh case is purely the
// multiply-scattered surface emission, so it directly tests the diffuse
// thermal field.
//-----------------------------------------------------------------------

void check_lidort_thermal_downwelling(boost::shared_ptr<LidortRtDriver>& ldriver,
                                      LidortDriverLambertianThermalFixture& t)
{
  BOOST_CHECK(ldriver->downwelling());
  double bb = planck(568.69, 290.0);
  Array<double, 1> bb_atm(t.nlayer + 1);
  bb_atm = bb;

  // (albedo, taur, taug, I_expected)
  double scen[3][4] = {
    {0.1, 1.0e-6, 1.0e-6, 2.0101793968e-07},
    {0.1, 1.0e-6, 1.0,    8.7632813980e-02},
    {0.9, 2.0e-2, 1.0e-6, 1.3981092762e-04} };
  for(int i = 0; i < 3; ++i) {
    t.surface_params(0) = scen[i][0];
    t.od = (scen[i][1] + scen[i][2]) / t.nlayer;
    t.ssa.value() = scen[i][1] / (scen[i][1] + scen[i][2]);
    Array<double, 1> refl(ldriver->reflectance_calculate(
        t.heights, t.sza(0), t.zen(0), t.azm(0), t.surface_type,
        t.surface_params, t.od.value(), t.ssa.value(), t.pf.value(), bb, bb_atm));
    BOOST_CHECK_CLOSE(refl(0), scen[i][3], 0.5);
  }
  // Gas case against the closed form as well
  BOOST_CHECK_CLOSE(scen[1][3], bb * (1 - exp(-1.0 - 1.0e-6)), 1e-4);
}

BOOST_AUTO_TEST_CASE(thermal_downwelling)
{
  Fixture1 t;
  BOOST_CHECK(!t.lidort_driver->downwelling());
  t.lidort_driver->set_downwelling(true);
  check_lidort_thermal_downwelling(t.lidort_driver, t);

  // The direction is carried in the serialized RT interface
  if(have_serialize_supported()) {
    boost::shared_ptr<LidortRtDriver> ldriver_r =
      serialize_read_string<LidortRtDriver>(serialize_write_string(t.lidort_driver));
    check_lidort_thermal_downwelling(ldriver_r, t);
  }

  // Switching back restores the upwelling (TOA) result
  t.lidort_driver->set_downwelling(false);
  t.run_test(t.lidort_driver);
}

// Atmospheric/surface blackbody, optical-depth and surface Jacobians vs
// finite differences for up/downwelling, black/reflecting surface and
// absorbing/scattering atmosphere. Regression test for two LIDORT LBBF bugs:
// the user-stream BRDF array was declared with transposed dimensions
// (missing most of the surface-reflected term in upwelling atmospheric BB
// Jacobians), and the 2-layer BVP used the wrong (upper) linear coefficient
// for the BOA level.
BOOST_AUTO_TEST_CASE(thermal_jacobian_matrix_fd)
{
  Fixture1 t;
  check_thermal_jacobians_fd(*t.lidort_driver, t.pf.value(), t.surface_type,
                             lambertian_test_surfaces(), {0});
  // nlayers > 2 takes the general-case BVP setup paths
  check_thermal_jacobians_fd(*t.lidort_driver, t.pf.value(), t.surface_type,
                             lambertian_test_surfaces(), {0}, true, true,
                             1e-2, 1e-6, 1e-6, 5);
  // Anisotropic scattering exercises the higher moments in the thermal
  // particular solution (VLIDORT was passed the wrong Pi matrix there)
  check_thermal_jacobians_fd(*t.lidort_driver, hg_phase_function(t.pf.value()), t.surface_type,
                             lambertian_test_surfaces(), {0});
}

// Thermal emission over a Cox-Munk surface: the emissivity (and its
// Jacobian) comes from Kirchhoff's law applied to the BRDF in the BRDF
// supplement, rather than 1 - albedo.
BOOST_AUTO_TEST_CASE(thermal_coxmunk_jacobian_fd)
{
  Fixture1 t;
  Array<double, 1> zen(1);
  zen = 0.0;
  LidortRtDriver d(t.nstreams, t.nmoms, false, COXMUNK, zen, false, false, true);
  d.lidort_interface()->lidort_modin().mbool().ts_do_deltam_scaling(false);
  d.set_plane_parallel();
  check_thermal_jacobians_fd(d, t.pf.value(), COXMUNK, coxmunk_test_surfaces(), {0, 1, 2});
}

BOOST_AUTO_TEST_SUITE_END()

///////////////////////////////////////////////////////////////////////////////////////

BOOST_FIXTURE_TEST_SUITE(lidort_driver_coxmunk_simple, GlobalFixture)

class Fixture1 : public LidortDriverCoxmunkFixture
{
public:
void run_test(boost::shared_ptr<LidortRtDriver>& ldriver) {
  Array<double, 1> refl_calc(1);
  Array<double, 1> refl_expected(1);

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
  BOOST_CHECK_MATRIX_CLOSE_TOL(refl_expected, refl_calc, 1e-3);

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

  Range ra = Range::all();

  Array<double, 1> refl_calc(1);
  Array<double, 1> refl_expt(1);

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

  blitz::Array<double, 3> jac_atm;
  blitz::Array<double, 2> jac_surf_param;
  blitz::Array<double, 1> jac_surf_temp;
  blitz::Array<double, 2> jac_atm_temp;

  ArrayAd<double, 1> lidort_surface(surface_params.shape(), 1);

  lidort_surface.value() = surface_params;
  lidort_surface.jacobian() = 1.0;

  ldriver->reflectance_and_jacobian_calculate(heights, sza(0), zen(0), azm(0),
                                                    surface_type, lidort_surface,
                                                    od, ssa, pf, refl_calc, 
                                                    jac_atm, jac_surf_param, jac_surf_temp, jac_atm_temp);

  BOOST_CHECK_EQUAL(check_brdf_inputs(lidort_driver), true);

  // Adjust analytic jacobians have same meaning as fd jacobians
  jac_surf_param(0, 0) *= lidort_surface.jacobian()(0,0);
  jac_surf_param(1, 0) *= lidort_surface.jacobian()(1,0);

  // Value for VLIDORT
  refl_expt = 0.70235315460259928;
  BOOST_CHECK_MATRIX_CLOSE_TOL(refl_expt, refl_calc, 1e-3);

  // Check surface jacobians against FD

  blitz::Array<double, 1> pert_values(surface_params.extent(firstDim)-1);
  pert_values = 1e-8, 1e-8, 1e-6;

  blitz::Array<double, 1> jac_surf_param_fd( jac_surf_param.rows() );
  Array<double, 1> refl_fd(1);

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

    jac_surf_param_fd(p_idx) = (refl_fd(0) - refl_calc(0)) / pert_values(p_idx);
  }

  BOOST_CHECK_MATRIX_CLOSE_TOL(jac_surf_param(ra, 0), jac_surf_param_fd, 1e-7);

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
  jac_surf_param(0, 0) *= lidort_surface.jacobian()(0,0);
  jac_surf_param(1, 0) *= lidort_surface.jacobian()(1,0);

  for(int p_idx = 0; p_idx < pert_values.extent(firstDim); p_idx++) {
    blitz::Array<double,1> surface_params_pert( surface_params.extent() );
    surface_params_pert = surface_params;
    surface_params_pert(p_idx) += pert_values(p_idx);

    refl_fd = ldriver->reflectance_calculate(heights, sza(0), zen(0), azm(0),
                                             surface_type, surface_params_pert,
                                             od.value(), ssa.value(), pf.value());

    jac_surf_param_fd(p_idx) = (refl_fd(0) - refl_calc(0)) / pert_values(p_idx);
  }

  BOOST_CHECK_MATRIX_CLOSE_TOL(jac_surf_param(ra, 0), jac_surf_param_fd, 1e-7);

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
  jac_surf_param(0, 0) *= lidort_surface.jacobian()(0,0);
  jac_surf_param(1, 0) *= lidort_surface.jacobian()(1,0);

  for(int p_idx = 0; p_idx < pert_values.extent(firstDim); p_idx++) {
    blitz::Array<double,1> surface_params_pert( surface_params.extent() );
    surface_params_pert = surface_params;
    surface_params_pert(p_idx) += pert_values(p_idx);

    refl_fd = ldriver->reflectance_calculate(heights, sza(0), zen(0), azm(0),
                                             surface_type, surface_params_pert,
                                             od.value(), ssa.value(), pf.value());

    jac_surf_param_fd(p_idx) = (refl_fd(0) - refl_calc(0)) / pert_values(p_idx);
  }

  BOOST_CHECK_MATRIX_CLOSE_TOL(jac_surf_param(ra, 0), jac_surf_param_fd, 1e-7);

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

  Range ra = Range::all();

  Array<double, 1> refl_calc(1);
  blitz::Array<double, 3> jac_atm;
  blitz::Array<double, 2> jac_surf_param;
  blitz::Array<double, 1> jac_surf_temp;
  blitz::Array<double, 2> jac_atm_temp;

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
  Array<double, 1> refl_expected(1); 
  refl_expected = 0.035435854422713485;
  BOOST_CHECK_MATRIX_CLOSE_TOL(refl_expected, refl_calc, 1e-3);

  // Check surface jacobians against FD

  blitz::Array<double, 1> pert_values(surface_params.rows());
  pert_values = 1e-8, 1e-8, 1e-6, 1e-6, 1e-6;

  blitz::Array<double, 1> jac_surf_param_fd( jac_surf_param.rows() );
  Array<double, 1> refl_fd(1);

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

    jac_surf_param_fd(p_idx) = (refl_fd(0) - refl_calc(0)) / pert_values(p_idx);
  }

  BOOST_CHECK_MATRIX_CLOSE_TOL(jac_surf_param(ra, 0), jac_surf_param_fd, 2e-7);

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

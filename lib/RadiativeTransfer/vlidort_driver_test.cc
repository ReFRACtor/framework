#include "vlidort_rt_driver.h"
#include "unit_test_support.h"
#include "vlidort_fixture.h"
#include "lidort_rt_driver.h"
#include "lidort_fixture.h"
#include "old_constant.h"
#include "fe_disable_exception.h"
#include "planck.h"

#include "spurr_brdf_types.h"

using namespace FullPhysics;
using namespace blitz;

bool check_brdf_inputs(boost::shared_ptr<VLidortRtDriver>& vlidort_driver) {
  VLidort_Vbrdf_Sup_Accessories vbrdf_check = VLidort_Vbrdf_Sup_Accessories(
      vlidort_driver->vlidort_interface()->vlidort_fixin_ptr(),
      vlidort_driver->vlidort_interface()->vlidort_modin_ptr(),
      vlidort_driver->vlidort_brdf_interface()->vbrdf_sup_in_ptr());
  vbrdf_check.vbrdf_input_check();

  VLidort_Exception_Handling& vbrdf_check_status = vbrdf_check.vlidort_vbrdfcheck_status();
  VLidort_Pars vlid_pars = VLidort_Pars::instance();

  if (vbrdf_check_status.ts_status_inputcheck() != vlid_pars.vlidort_success())
    std::cerr << vbrdf_check_status << std::endl;

  return vbrdf_check_status.ts_status_inputcheck() == vlid_pars.vlidort_success();
}

//=======================================================================
// nstokes = 1 (scalar). The intensity should match the values computed
// by LIDORT for the equivalent scalar problem.
//=======================================================================

BOOST_FIXTURE_TEST_SUITE(vlidort_driver_lambertian_solar, GlobalFixture)

class Fixture1 : public VLidortDriverLambertianFixture
{
public:
  Fixture1() : VLidortDriverLambertianFixture(1) {}

  void run_test(boost::shared_ptr<VLidortRtDriver>& driver) {
    Array<double, 1> refl_calc;
    Array<double, 1> refl_expected(1);

    ////////////////
    // Surface only
    surface_params(0) = 1.0;

    taur = 1.0e-6/nlayer;
    taug = 1.0e-6/nlayer;

    od = taur + taug;
    ssa.value() = taur / od.value();

    refl_calc.reference(driver->reflectance_calculate(heights, sza(0), zen(0), azm(0),
                                                      surface_type, surface_params,
                                                      od.value(), ssa.value(), pf.value()));

    // Surface only = 1/pi. Tolerance is a percentage difference.
    refl_expected = 1/OldConstant::pi * surface_params(0);
    BOOST_CHECK_MATRIX_CLOSE_TOL(refl_expected, refl_calc, 1e-3);

    ////////////////
    // Gas + Surface
    surface_params(0) = 1.0;

    taur = 1.0e-6/nlayer;
    taug = 1.0/nlayer;

    od = taur + taug;
    ssa.value() = taur / od.value();

    refl_calc.reference(driver->reflectance_calculate(heights, sza(0), zen(0), azm(0),
                                                      surface_type, surface_params,
                                                      od.value(), ssa.value(), pf.value()));

    refl_expected = 1/OldConstant::pi * exp(-1/cos(sza(0))) * exp(-1/cos(zen(0)));
    BOOST_CHECK_MATRIX_CLOSE_TOL(refl_expected, refl_calc, 1e-3);

    ////////////////
    // Rayleigh only
    surface_params(0) = 0.0;

    taur = 2.0e-2/nlayer;
    taug = 1.0e-6/nlayer;

    od = taur + taug;
    ssa.value() = taur / od.value();

    refl_calc.reference(driver->reflectance_calculate(heights, sza(0), zen(0), azm(0),
                                                      surface_type, surface_params,
                                                      od.value(), ssa.value(), pf.value()));

    // Expected value from VLIDORT, this is not going to agree to infinite
    // precision. Tolerance is a percentage.
    refl_expected = 2.387246757232095E-003;
    BOOST_CHECK_MATRIX_CLOSE_TOL(refl_expected, refl_calc, 6e-3);
  }
};

BOOST_AUTO_TEST_CASE(solar_sources)
{
  Fixture1 t;
  t.run_test(t.vlidort_driver);
}

BOOST_AUTO_TEST_CASE(serialization)
{
  if(!have_serialize_supported())
    return;
  Fixture1 t;
  std::string d = serialize_write_string(t.vlidort_driver);
  if(false)
    std::cerr << d;
  boost::shared_ptr<VLidortRtDriver> driver_r =
    serialize_read_string<VLidortRtDriver>(d);
  t.run_test(driver_r);
}

BOOST_AUTO_TEST_SUITE_END()

//=======================================================================
// Thermal emission only (no solar source), nstokes = 1. The Stokes-I
// radiance is checked against the offline VLIDORT thermal reference
// (reference/rt/vlidort_thermal_test, results.1.log) for the same three
// scenarios used by the scalar LIDORT thermal test. The blackbody
// weighting-function Jacobians (jac_surf_temp, jac_atm_temp) are checked
// against finite differences.
//=======================================================================

BOOST_FIXTURE_TEST_SUITE(vlidort_driver_lambertian_thermal, GlobalFixture)

class Fixture1 : public VLidortDriverLambertianThermalFixture
{
public:
  Fixture1() : VLidortDriverLambertianThermalFixture(1) {}

  void run_test(boost::shared_ptr<VLidortRtDriver>& driver) {
    Array<double, 1> refl_calc(1);
    Array<double, 1> refl_expt(1);

    /////////////////////////////////////////////
    // Thermal setup

    // Just some nominal values to calculate the planck function. These
    // match the scalar LIDORT thermal test and the offline VLIDORT
    // reference driver.
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

    refl_calc = driver->reflectance_calculate(heights, sza(0), azm(0), zen(0),
                                              surface_type, surface_params,
                                              od.value(), ssa.value(), pf.value(),
                                              bb_surface, bb_atm);

    // Expected values calculated offline in VLIDORT (results.1.log)
    refl_expt = 0.12476970639470732;
    BOOST_CHECK_MATRIX_CLOSE_TOL(refl_expt, refl_calc, 1e-3);

    ////////////////////////
    // Thermal gas + surface
    surface_params(0) = 0.1;

    taur = 1.0e-6/nlayer;
    taug = 1.0/nlayer;

    od = taur + taug;
    ssa.value() = taur / od.value();

    refl_calc = driver->reflectance_calculate(heights, sza(0), azm(0), zen(0),
                                              surface_type, surface_params,
                                              od.value(), ssa.value(), pf.value(),
                                              bb_surface, bb_atm);

    refl_expt = 0.13751414471789294;
    BOOST_CHECK_MATRIX_CLOSE_TOL(refl_expt, refl_calc, 1e-3);

    /////////////////////////////////////////////
    // Thermal rayleigh only, no surface, no gas
    surface_params(0) = 0.9;

    taur = 2.0e-2/nlayer;
    taug = 1.0e-6/nlayer;

    od = taur + taug;
    ssa.value() = taur / od.value();

    refl_calc = driver->reflectance_calculate(heights, sza(0), azm(0), zen(0),
                                              surface_type, surface_params,
                                              od.value(), ssa.value(), pf.value(),
                                              bb_surface, bb_atm);

    refl_expt = 1.3966249402729437e-2;
    BOOST_CHECK_MATRIX_CLOSE_TOL(refl_expt, refl_calc, 7e-3);
  }
};

BOOST_AUTO_TEST_CASE(thermal_emission)
{
  Fixture1 t;
  t.run_test(t.vlidort_driver);
}

// Blackbody weighting-function Jacobians (jac_surf_temp = d(stokes)/d(surface
// blackbody), jac_atm_temp = d(stokes)/d(per-level atmosphere blackbody))
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

  t.vlidort_driver->reflectance_and_jacobian_calculate(
      t.heights, t.sza(0), t.azm(0), t.zen(0), t.surface_type, surface,
      od, ssa, pf, refl, jac_atm, jac_surf_param, jac_surf_temp, jac_atm_temp,
      bb_surface, bb_atm);

  double pert = 1e-3;

  // Surface blackbody jacobian
  {
    Array<double, 1> refl_pert(t.vlidort_driver->reflectance_calculate(
        t.heights, t.sza(0), t.azm(0), t.zen(0), t.surface_type, surface.value(),
        od.value(), ssa.value(), pf.value(), bb_surface + pert, bb_atm));
    for(int s = 0; s < refl.rows(); ++s) {
      double jac_fd = (refl_pert(s) - refl(s)) / pert;
      BOOST_CHECK_CLOSE(jac_surf_temp(s), jac_fd, 1e-2);
    }
  }

  // Per-level atmosphere blackbody jacobian
  for(int lev = 0; lev < bb_atm.rows(); ++lev) {
    Array<double, 1> bb_atm_pert(bb_atm.copy());
    bb_atm_pert(lev) += pert;
    Array<double, 1> refl_pert(t.vlidort_driver->reflectance_calculate(
        t.heights, t.sza(0), t.azm(0), t.zen(0), t.surface_type, surface.value(),
        od.value(), ssa.value(), pf.value(), bb_surface, bb_atm_pert));
    for(int s = 0; s < refl.rows(); ++s) {
      double jac_fd = (refl_pert(s) - refl(s)) / pert;
      BOOST_CHECK_CLOSE(jac_atm_temp(lev, s), jac_fd, 1e-2);
    }
  }
}

// Cross-solver consistency: for the scalar (nstokes = 1) thermal problem
// VLIDORT and LIDORT wrap equivalent discrete-ordinate physics, so the VLIDORT
// Stokes-I radiance must match the LIDORT reflectance for the same
// configuration. Both fixtures share geometry, streams, layering and surface.
BOOST_AUTO_TEST_CASE(lidort_consistency)
{
  Fixture1 vt;                                 // VLIDORT nstokes=1 thermal
  LidortDriverLambertianThermalFixture lt;     // LIDORT thermal, same config

  double wn = 568.69;
  double temperature = 290.0;
  double bb_surface = planck(wn, temperature);
  Array<double, 1> bb_atm(vt.nlayer + 1);
  bb_atm = planck(wn, temperature);

  // (surface albedo, taur, taug) for the three thermal scenarios.
  double scen[3][3] = { {0.1, 1.0e-6, 1.0e-6},
                        {0.1, 1.0e-6, 1.0},
                        {0.9, 2.0e-2, 1.0e-6} };

  for(int i = 0; i < 3; ++i) {
    double albedo = scen[i][0];
    double taur = scen[i][1]/vt.nlayer;
    double taug = scen[i][2]/vt.nlayer;
    double od_v = taur + taug;
    double ssa_v = taur / od_v;

    Array<double, 1> od1(1), ssa1(1);
    od1 = od_v; ssa1 = ssa_v;

    vt.surface_params(0) = albedo;
    Array<double, 1> refl_v(vt.vlidort_driver->reflectance_calculate(
        vt.heights, vt.sza(0), vt.azm(0), vt.zen(0), vt.surface_type,
        vt.surface_params, od1, ssa1, vt.pf.value(), bb_surface, bb_atm));

    lt.surface_params(0) = albedo;
    Array<double, 1> refl_l(lt.lidort_driver->reflectance_calculate(
        lt.heights, lt.sza(0), lt.zen(0), lt.azm(0), lt.surface_type,
        lt.surface_params, od1, ssa1, lt.pf.value(), bb_surface, bb_atm));

    BOOST_CHECK_MATRIX_CLOSE_TOL(refl_l, refl_v, 1e-7);
  }
}

BOOST_AUTO_TEST_CASE(serialization)
{
  if(!have_serialize_supported())
    return;
  Fixture1 t;
  std::string d = serialize_write_string(t.vlidort_driver);
  if(false)
    std::cerr << d;
  boost::shared_ptr<VLidortRtDriver> driver_r =
    serialize_read_string<VLidortRtDriver>(d);
  t.run_test(driver_r);
}

BOOST_AUTO_TEST_SUITE_END()

//=======================================================================
// nstokes = 1, 3 and 4 (vector) in an off-nadir polarizing geometry. The
// full Stokes vector is checked against the offline VLIDORT reference
// driver (test_vlidort_refractor.f90 in
// reference/rt/vlidort_2p8p3_lin_tester), for each sphericity mode
// (plane-parallel, pseudo-spherical, line-of-sight). Pure Rayleigh over a
// black surface so the comparison isolates the vector scattering /
// greek-matrix mapping.
//=======================================================================

BOOST_FIXTURE_TEST_SUITE(vlidort_driver_vector, GlobalFixture)

enum Sphericity { PLANE_PARALLEL, PSEUDO_SPHERICAL, LINE_OF_SIGHT };

void set_sphericity(boost::shared_ptr<VLidortRtDriver>& driver, Sphericity mode) {
  switch(mode) {
  case PLANE_PARALLEL:   driver->set_plane_parallel(); break;
  case PSEUDO_SPHERICAL: driver->set_pseudo_spherical(); break;
  case LINE_OF_SIGHT:    driver->set_line_of_sight(); break;
  }
}

void run_vector_test(VLidortDriverVectorFixture& t, Sphericity mode,
                     const blitz::Array<double, 1>& stokes_expected, double tol) {
  set_sphericity(t.vlidort_driver, mode);

  // Pure Rayleigh atmosphere, black surface (albedo 0).
  t.surface_params(0) = 0.0;

  t.taur = 2.0e-2/t.nlayer;
  t.taug = 1.0e-6/t.nlayer;

  t.od = t.taur + t.taug;
  t.ssa.value() = t.taur / t.od.value();

  // reflectance_calculate signature is (height, sza, azm, zen, ...)
  Array<double, 1> stokes_calc(t.vlidort_driver->reflectance_calculate(
      t.heights, t.sza(0), t.azm(0), t.zen(0),
      t.surface_type, t.surface_params,
      t.od.value(), t.ssa.value(), t.pf.value()));

  BOOST_CHECK_EQUAL(stokes_calc.rows(), stokes_expected.rows());
  BOOST_CHECK_MATRIX_CLOSE_TOL(stokes_expected, stokes_calc, tol);
}

// ---- Plane-parallel (reference mode 0) ----

BOOST_AUTO_TEST_CASE(plane_parallel_nstokes_1)
{
  VLidortDriverVectorFixture t(1);
  Array<double, 1> e(1);
  e = 1.6883136744104341e-03;
  run_vector_test(t, PLANE_PARALLEL, e, 1e-7);
}

BOOST_AUTO_TEST_CASE(plane_parallel_nstokes_3)
{
  VLidortDriverVectorFixture t(3);
  Array<double, 1> e(3);
  e = 1.6821888550733900e-03, -1.5663936555950202e-04, -1.2398207619582417e-03;
  run_vector_test(t, PLANE_PARALLEL, e, 1e-7);
}

BOOST_AUTO_TEST_CASE(plane_parallel_nstokes_4)
{
  VLidortDriverVectorFixture t(4);
  Array<double, 1> e(4);
  e = 1.6821888550731530e-03, -1.5663936556101861e-04, -1.2398207619582389e-03, 0.0;
  run_vector_test(t, PLANE_PARALLEL, e, 1e-7);
}

// ---- Pseudo-spherical (reference mode 1) ----

BOOST_AUTO_TEST_CASE(pseudo_spherical_nstokes_1)
{
  VLidortDriverVectorFixture t(1);
  Array<double, 1> e(1);
  e = 1.6884936368769267e-03;
  run_vector_test(t, PSEUDO_SPHERICAL, e, 1e-7);
}

BOOST_AUTO_TEST_CASE(pseudo_spherical_nstokes_3)
{
  VLidortDriverVectorFixture t(3);
  Array<double, 1> e(3);
  e = 1.6823681621588628e-03, -1.5665606058624514e-04, -1.2399529087074581e-03;
  run_vector_test(t, PSEUDO_SPHERICAL, e, 1e-7);
}

// ---- Line-of-sight / outgoing (reference mode 2) ----

BOOST_AUTO_TEST_CASE(line_of_sight_nstokes_1)
{
  VLidortDriverVectorFixture t(1);
  Array<double, 1> e(1);
  e = 1.6823529986010875e-03;
  run_vector_test(t, LINE_OF_SIGHT, e, 1e-7);
}

BOOST_AUTO_TEST_CASE(line_of_sight_nstokes_3)
{
  VLidortDriverVectorFixture t(3);
  Array<double, 1> e(3);
  e = 1.6762275238830235e-03, -1.5607360168268711e-04, -1.2353734844838193e-03;
  run_vector_test(t, LINE_OF_SIGHT, e, 1e-7);
}

BOOST_AUTO_TEST_SUITE_END()

//=======================================================================
// Analytic jacobians (the linearized greek / copy_jacobians path) checked
// against finite differences of the intensity for the full Stokes vector.
//=======================================================================

BOOST_FIXTURE_TEST_SUITE(vlidort_driver_vector_jacobian, GlobalFixture)

// Atmospheric (optical depth) profile jacobian versus finite difference.
BOOST_AUTO_TEST_CASE(atmospheric_jacobian_fd)
{
  FeDisableException disable_fp;
  VLidortDriverVectorFixture t(3);
  t.vlidort_driver->set_plane_parallel();

  double taur = 2.0e-2, taug = 1.0e-6;

  // od carries a single derivative d(od)/dx = 1; ssa and pf held fixed so
  // the analytic jacobian is d(stokes)/d(od).
  ArrayAd<double, 1> od(1, 1), ssa(1, 1);
  od.value()(0) = taur + taug;
  od.jacobian()(0, 0) = 1.0;
  ssa.value()(0) = taur / od.value()(0);
  ssa.jacobian()(0, 0) = 0.0;

  ArrayAd<double, 3> pf(t.pf);
  pf.jacobian() = 0.0;

  ArrayAd<double, 1> surface(1, 1);
  surface.value()(0) = 0.0;   // black surface
  surface.jacobian()(0, 0) = 0.0;

  Array<double, 1> refl;
  Array<double, 3> jac_atm;
  Array<double, 2> jac_surf_param, jac_atm_temp;
  Array<double, 1> jac_surf_temp;

  t.vlidort_driver->reflectance_and_jacobian_calculate(
      t.heights, t.sza(0), t.azm(0), t.zen(0), t.surface_type, surface,
      od, ssa, pf, refl, jac_atm, jac_surf_param, jac_surf_temp, jac_atm_temp);

  // Finite difference: perturb od value only.
  double pert = 1e-7;
  Array<double, 1> od_pert(od.value().copy());
  od_pert(0) += pert;
  Array<double, 1> refl_pert(t.vlidort_driver->reflectance_calculate(
      t.heights, t.sza(0), t.azm(0), t.zen(0), t.surface_type, surface.value(),
      od_pert, ssa.value(), pf.value()));

  for(int s = 0; s < refl.rows(); ++s) {
    double jac_fd = (refl_pert(s) - refl(s)) / pert;
    // Absolute tolerance: the tiny Q/U components sit near finite-difference
    // roundoff, where a relative check is meaningless.
    BOOST_CHECK_SMALL(jac_atm(0, 0, s) - jac_fd, 1e-5);
  }
}

// Surface (Lambertian albedo) jacobian versus finite difference.
BOOST_AUTO_TEST_CASE(surface_jacobian_fd)
{
  FeDisableException disable_fp;
  VLidortDriverVectorFixture t(3);
  t.vlidort_driver->set_plane_parallel();

  double taur = 2.0e-2, taug = 1.0e-6;
  t.od = taur + taug;
  t.ssa.value() = taur / t.od.value();

  ArrayAd<double, 1> surface(1, 1);
  surface.value()(0) = 0.1;
  surface.jacobian()(0, 0) = 1.0;

  Array<double, 1> refl;
  Array<double, 3> jac_atm;
  Array<double, 2> jac_surf_param, jac_atm_temp;
  Array<double, 1> jac_surf_temp;

  t.vlidort_driver->reflectance_and_jacobian_calculate(
      t.heights, t.sza(0), t.azm(0), t.zen(0), t.surface_type, surface,
      t.od, t.ssa, t.pf, refl, jac_atm, jac_surf_param, jac_surf_temp, jac_atm_temp);

  // Finite difference: perturb the albedo.
  double pert = 1e-7;
  Array<double, 1> surface_pert(surface.value().copy());
  surface_pert(0) += pert;
  Array<double, 1> refl_pert(t.vlidort_driver->reflectance_calculate(
      t.heights, t.sza(0), t.azm(0), t.zen(0), t.surface_type, surface_pert,
      t.od.value(), t.ssa.value(), t.pf.value()));

  for(int s = 0; s < refl.rows(); ++s) {
    double jac_fd = (refl_pert(s) - refl(s)) / pert;
    BOOST_CHECK_SMALL(jac_surf_param(0, s) - jac_fd, 1e-5);
  }
}

BOOST_AUTO_TEST_SUITE_END()

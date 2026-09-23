#include "vlidort_rt_driver.h"
#include "unit_test_support.h"
#include "vlidort_fixture.h"
#include "old_constant.h"

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

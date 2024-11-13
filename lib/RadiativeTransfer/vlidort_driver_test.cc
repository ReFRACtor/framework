#include "vlidort_rt_driver.h"
#include "unit_test_support.h"
#include "lidort_fixture.h"
#include "old_constant.h"

#include "spurr_brdf_types.h"

using namespace FullPhysics;
using namespace blitz;

bool check_brdf_inputs(boost::shared_ptr<VLidortRtDriver>& vlidort_driver) {
  VLidort_VBrdf_Sup_Accessories vbrdf_check = VLidort_VBrdf_Sup_Accessories(
      vlidort_driver->vlidort_interface()->vlidort_fixin_ptr(),
      vlidort_driver->vlidort_interface()->vlidort_modin_ptr(),
      vlidort_driver->vlidort_brdf_interface()->vbrdf_sup_in_ptr());
  vbrdf_check.vbrdf_input_check();

  VLidort_Exception_Handling& vbrdf_check_status = vbrdf_check.vlidort_brdfcheck_status();
  VLidort_Pars lid_pars = VLidort_Pars::instance();

  if (vbrdf_check_status.ts_status_inputcheck() != lid_pars.vlidort_success())
    std::cerr << vbrdf_check_status << std::endl;

  return vbrdf_check_status.ts_status_inputcheck() == vlid_pars.lidort_success();
}

BOOST_FIXTURE_TEST_SUITE(vlidort_driver_lambertian_solar, GlobalFixture)

class Fixture1 : public VLidortDriverLambertianFixture
{
public:
void run_test(boost::shared_ptr<VLidortRtDriver>& driver) {
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

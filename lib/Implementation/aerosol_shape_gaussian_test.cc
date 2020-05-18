#include "aerosol_shape_gaussian.h"
#include "met_data_fixture.h"
#include "unit_test_support.h"

using namespace FullPhysics;
using namespace blitz;

BOOST_FIXTURE_TEST_SUITE(aerosol_shape_gaussian, MetDataFixture)

BOOST_AUTO_TEST_CASE(basic)
{
  // Arrays defining coeffs
  blitz::Array<bool, 1> ret_flag(3);
  ret_flag = true;
  blitz::Array<double, 1> coeffs(3);

  // Loaded expected results
  IfstreamCs gaussian_expt(test_data_dir() + "expected/aerosol_shape/gaussian");

  // Kahn
  Array<double, 1> kahn_expt;
  gaussian_expt >> kahn_expt;

  coeffs = -3.28341, 1.0, 0.2;
  AerosolShapeGaussian aer_kahn_log = AerosolShapeGaussian(pressure, ret_flag, coeffs, "Kahn", false);
  BOOST_CHECK_MATRIX_CLOSE_TOL(kahn_expt, aer_kahn_log.aerosol_extinction().value(), 1e-14);

  coeffs(0) = exp(coeffs(0));
  AerosolShapeGaussian aer_kahn_lin = AerosolShapeGaussian(pressure, ret_flag, coeffs, "Kahn", true);
  BOOST_CHECK_MATRIX_CLOSE_TOL(kahn_expt, aer_kahn_lin.aerosol_extinction().value(), 1e-14);

  // Water
  Array<double, 1> water_expt;
  gaussian_expt >> water_expt;

  coeffs = -3.28341, 0.75, 0.1;
  AerosolShapeGaussian aer_water_log = AerosolShapeGaussian(pressure, ret_flag, coeffs, "Water", false);
  BOOST_CHECK_MATRIX_CLOSE_TOL(water_expt, aer_water_log.aerosol_extinction().value(), 1e-14);

  coeffs(0) = exp(coeffs(0));
  AerosolShapeGaussian aer_water_lin = AerosolShapeGaussian(pressure, ret_flag, coeffs, "Water", true);
  BOOST_CHECK_MATRIX_CLOSE_TOL(water_expt, aer_water_lin.aerosol_extinction().value(), 1e-14);

  // Ice
  Array<double, 1> ice_expt;
  gaussian_expt >> ice_expt;

  coeffs = -3.28341, 0.3, 0.04;
  AerosolShapeGaussian aer_ice_log = AerosolShapeGaussian(pressure, ret_flag, coeffs, "Ice", false);
  BOOST_CHECK_MATRIX_CLOSE_TOL(ice_expt, aer_ice_log.aerosol_extinction().value(), 1e-14);

  coeffs(0) = exp(coeffs(0));
  AerosolShapeGaussian aer_ice_lin = AerosolShapeGaussian(pressure, ret_flag, coeffs, "Ice", true);
  BOOST_CHECK_MATRIX_CLOSE_TOL(ice_expt, aer_ice_lin.aerosol_extinction().value(), 1e-14);
  
}

BOOST_AUTO_TEST_SUITE_END()


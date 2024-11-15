#include "aerosol_extinction_linear.h"
#include "met_data_fixture.h"
#include "unit_test_support.h"

using namespace FullPhysics;
using namespace blitz;

BOOST_FIXTURE_TEST_SUITE(aerosol_extinction_linear, MetDataFixture)

BOOST_AUTO_TEST_CASE(basic)
{
  // Arrays defining coeffs
  blitz::Array<double, 1> coeffs(3);

  // Kahn
  Array<double, 1> kahn_linear(3);
  coeffs = 0.0, 1.0, 0.2;
  kahn_linear = coeffs;

  AerosolExtinctionLinear aer_kahn_lin = AerosolExtinctionLinear(pressure, coeffs, "Kahn");
  BOOST_CHECK_MATRIX_CLOSE_TOL(kahn_linear, aer_kahn_lin.aerosol_extinction().value(), 1e-14);
  BOOST_CHECK_MATRIX_CLOSE_TOL(kahn_linear, aer_kahn_lin.aerosol_extinction(Pressure::NATIVE_ORDER).value(), 1e-14);
  BOOST_CHECK_MATRIX_CLOSE_TOL(kahn_linear.reverse(blitz::firstDim), aer_kahn_lin.aerosol_extinction(Pressure::DECREASING_PRESSURE).value(), 1e-14);

  // Water
  Array<double, 1> water_expt(3);
  coeffs = 10.0, 0.75, 0.1;
  water_expt = coeffs;

  AerosolExtinctionLinear aer_water_lin = AerosolExtinctionLinear(pressure_reverse, coeffs.reverse(blitz::firstDim), "Water");
  BOOST_CHECK_MATRIX_CLOSE_TOL(water_expt, aer_water_lin.aerosol_extinction().value(), 1e-14);
  BOOST_CHECK_MATRIX_CLOSE_TOL(water_expt.reverse(blitz::firstDim), aer_water_lin.aerosol_extinction(Pressure::NATIVE_ORDER).value(), 1e-14);
  BOOST_CHECK_MATRIX_CLOSE_TOL(water_expt.reverse(blitz::firstDim), aer_water_lin.aerosol_extinction(Pressure::DECREASING_PRESSURE).value(), 1e-14);

  // Ice
  Array<double, 1> ice_expt(3);
  coeffs = 8.0, 0.3, 0.04;
  ice_expt = coeffs;

  AerosolExtinctionLinear aer_ice_lin = AerosolExtinctionLinear(pressure, coeffs, "Ice");
  BOOST_CHECK_MATRIX_CLOSE_TOL(ice_expt, aer_ice_lin.aerosol_extinction().value(), 1e-14);
  
}

BOOST_AUTO_TEST_SUITE_END()


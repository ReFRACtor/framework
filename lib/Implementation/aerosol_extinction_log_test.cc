#include "aerosol_extinction_log.h"
#include "met_data_fixture.h"
#include "unit_test_support.h"

using namespace FullPhysics;
using namespace blitz;

BOOST_FIXTURE_TEST_SUITE(aerosol_extinction_log, MetDataFixture)

BOOST_AUTO_TEST_CASE(basic)
{
  // Arrays defining coeffs
  blitz::Array<bool, 1> ret_flag(3);
  ret_flag = true;
  blitz::Array<double, 1> coeffs(3);

  // Kahn
  Array<double, 1> kahn_log(3);
  coeffs = 0.0, 1.0, 0.2;
  kahn_log = exp(coeffs);

  AerosolExtinctionLog aer_kahn_log = AerosolExtinctionLog(pressure, ret_flag, coeffs, "Kahn");
  BOOST_CHECK_MATRIX_CLOSE_TOL(kahn_log, aer_kahn_log.aerosol_extinction().value(), 1e-14);

  // Water
  Array<double, 1> water_expt(3);
  coeffs = 10.0, 0.75, 0.1;
  water_expt = exp(coeffs);

  AerosolExtinctionLog aer_water_log = AerosolExtinctionLog(pressure, ret_flag, coeffs, "Water");
  BOOST_CHECK_MATRIX_CLOSE_TOL(water_expt, aer_water_log.aerosol_extinction().value(), 1e-14);

  // Ice
  Array<double, 1> ice_expt(3);
  coeffs = 8.0, 0.3, 0.04;
  ice_expt = exp(coeffs);

  AerosolExtinctionLog aer_ice_log = AerosolExtinctionLog(pressure, ret_flag, coeffs, "Ice");
  BOOST_CHECK_MATRIX_CLOSE_TOL(ice_expt, aer_ice_log.aerosol_extinction().value(), 1e-14);
  
}

BOOST_AUTO_TEST_SUITE_END()


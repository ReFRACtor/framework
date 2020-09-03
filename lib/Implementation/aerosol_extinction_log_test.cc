#include "aerosol_extinction_log.h"
#include "met_data_fixture.h"
#include "unit_test_support.h"

using namespace FullPhysics;
using namespace blitz;

BOOST_FIXTURE_TEST_SUITE(aerosol_extinction_log, MetDataFixture)

BOOST_AUTO_TEST_CASE(basic)
{
  // Arrays defining coeffs
  blitz::Array<double, 1> aext(20);

  aext = 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20;

  // Check that the fm_view of aext matches what we put in
  AerosolExtinctionLog aer_water_log = AerosolExtinctionLog(pressure, aext, "Water");
  BOOST_CHECK_MATRIX_CLOSE_TOL(aer_water_log.aerosol_extinction().value(), aext, 1e-14);

  // Check that coefficients are stored in log within SubStateVectorArray
  for(int j = 0; j < aext.rows(); j++) {
      BOOST_CHECK_CLOSE(aer_water_log.coefficient().value()(j), log(aext(j)), 1e-8);
  }
  
}

BOOST_AUTO_TEST_SUITE_END()


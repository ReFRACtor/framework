#include "unit_test_support.h"
#include "atmosphere_fixture.h"
#include "relative_humidity.h"

using namespace FullPhysics;
using namespace blitz;
BOOST_FIXTURE_TEST_SUITE(relative_humidity, AtmosphereFixture)

BOOST_AUTO_TEST_CASE(basic)
{
  RelativeHumidity h(atm->absorber_ptr(), atm->temperature_ptr(),
                     atm->pressure_ptr());
  BOOST_CHECK_CLOSE(h.specific_humidity_grid()(18).value(), 0.0071139603422119601, 1e-2);
  BOOST_CHECK_CLOSE(h.relative_humidity_grid()(18).value(), 80.827494397125776, 1e-2);
}

BOOST_AUTO_TEST_SUITE_END()

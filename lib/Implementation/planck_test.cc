#include "unit_test_support.h"
#include "planck.h"
#include "auto_derivative.h"

using namespace FullPhysics;

BOOST_FIXTURE_TEST_SUITE(planck_function, GlobalFixture)

BOOST_AUTO_TEST_CASE(basic)
{
    // Test against value from:
    // https://ncc.nesdis.noaa.gov/data/planck.html
    // Note their units are in mW/m^2/sr/cm^-1, so divide by 1000 to get in W/...
    double expected_bb = 0.09953544145262894;

    double wn = 908.8080;
    double temperature = 290.0;
    AutoDerivative<double> bbody = planck(wn, temperature);

    BOOST_CHECK_CLOSE_FRACTION(bbody.value(), expected_bb, 1e-4);
}

BOOST_AUTO_TEST_SUITE_END()

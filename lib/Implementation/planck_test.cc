#include "unit_test_support.h"
#include "planck.h"
#include "auto_derivative.h"

using namespace FullPhysics;
using namespace blitz;

BOOST_FIXTURE_TEST_SUITE(planck_function, GlobalFixture)

BOOST_AUTO_TEST_CASE(basic)
{
    // Test against value from:
    // https://ncc.nesdis.noaa.gov/data/planck.html
    // Note their units are in mW/m^2/sr/cm^-1, so divide by 1000 to get in W/...
    double expected_bb = 0.09953544145262894;

    double wn = 908.8080;
    double temperature = 290.0;
    double bbody = planck(wn, temperature);

    BOOST_CHECK_CLOSE_FRACTION(bbody, expected_bb, 1e-4);
}

BOOST_AUTO_TEST_CASE(jacobian)
{
    double expected_bb = 0.09953544145262894;

    // Value from previous run, assume get_planck is well tested
    double expected_grad = 0.00156477;

    double wn = 908.8080;
    double temperature = 290.0;

    Array<double, 1> gradient(1);
    gradient = 1;

    AutoDerivative<double> bbody = planck(wn, temperature, gradient);

    BOOST_CHECK_CLOSE_FRACTION(bbody.value(), expected_bb, 1e-4);
    BOOST_CHECK_CLOSE_FRACTION(bbody.gradient()(0), expected_grad, 1e-4);
}

BOOST_AUTO_TEST_SUITE_END()

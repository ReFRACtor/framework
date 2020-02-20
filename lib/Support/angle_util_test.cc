#include "unit.h"
#include "unit_test_support.h"

#include "angle_util.h"

using namespace FullPhysics;

BOOST_FIXTURE_TEST_SUITE(angle_util, GlobalFixture)

BOOST_AUTO_TEST_CASE(basic)
{
    // Angles come out of test case where I was testing OMI against MUSES
    DoubleWithUnit vza(50.36930465698242, units::deg);
    DoubleWithUnit sza(30.606971740722656, units::deg);
    DoubleWithUnit raz(175.28192138671875, units::deg);

    DoubleWithUnit sca = scattering_angle(vza, sza, raz);

    BOOST_CHECK_CLOSE(sca.convert(units::deg).value, 19.9862759301103, 1e-10);
}

BOOST_AUTO_TEST_SUITE_END()


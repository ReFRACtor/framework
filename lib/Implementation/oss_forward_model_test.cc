#include "oss_forward_model.h"
#include "unit_test_support.h"
#include "atmosphere_fixture.h"

using namespace FullPhysics;

BOOST_FIXTURE_TEST_SUITE(oss_forward_constant, AtmosphereFixture)

BOOST_AUTO_TEST_CASE(radiance)
{
  OssForwardModel fm(atm);
  fm.setup_grid();
  // TODO: Add checks
}

BOOST_AUTO_TEST_SUITE_END()

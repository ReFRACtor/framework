#include "absorber_vmr_level.h"
#include "met_data_fixture.h"
#include "unit_test_support.h"

using namespace FullPhysics;
using namespace blitz;

BOOST_FIXTURE_TEST_SUITE(absorber_vmr_level, MetDataFixture)

BOOST_AUTO_TEST_CASE(basic)
{
  Array<double, 1> vmr(19);
  vmr = 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19;
  AbsorberVmrLevel avmr(pressure, vmr, "CO2");
  for(int i = 0; i < pressure->pressure_grid().rows(); ++i)
    BOOST_CHECK_CLOSE(avmr.volume_mixing_ratio(pressure->pressure_grid()(i).value).value(), vmr(i), 1e-8);
}

BOOST_AUTO_TEST_SUITE_END()

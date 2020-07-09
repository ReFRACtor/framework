#include "temperature_level.h"
#include "unit_test_support.h"
#include "met_data_fixture.h"

using namespace FullPhysics;
using namespace blitz;

BOOST_FIXTURE_TEST_SUITE(temperature_level, MetDataFixture)

BOOST_AUTO_TEST_CASE(basic)
{
  StateVector sv;

  Array<double, 1> temp_expect(19);
  temp_expect = 244.2, 214.553, 218.029, 222.544, 218.341, 221.37, 227.38,
    233.493, 239.376, 244.52, 248.708, 251.979, 254.537, 256.655, 258.521,
    260.155, 261.747, 261.732, 258.598;

  Array<bool, 1> temp_flag(temp_expect.rows());
  temp_flag = true;

  TemperatureLevel temp_lev(temp_expect, temp_flag, pressure);
  sv.add_observer(temp_lev);

  for(int i = 0; i < temp_expect.rows(); ++i) {
    BOOST_CHECK_CLOSE(temp_lev.temperature(pressure->pressure_grid()(i)).convert(units::K).value.value(), temp_expect(i), 1e-3);
  }
}

BOOST_AUTO_TEST_SUITE_END()


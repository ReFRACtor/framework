#include "temperature_met.h"
#include "met_data_fixture.h"
#include "pressure_sigma.h"
#include "unit_test_support.h"
#include "temperature_level_offset.h"

using namespace FullPhysics;
using namespace blitz;

BOOST_FIXTURE_TEST_SUITE(temperature_level_offset, MetDataFixture)

BOOST_AUTO_TEST_CASE(basic)
{
  boost::shared_ptr<Pressure> pressure_sigma(new PressureSigma(pressure_in, psurf_in, false));
  StateVector sv;

  Array<double, 1> temp_expect(19);
  temp_expect = 244.2, 214.553, 218.029, 222.544, 218.341, 221.37, 227.38,
    233.493, 239.376, 244.52, 248.708, 251.979, 254.537, 256.655, 258.521,
    260.155, 261.747, 261.732, 258.598;

  TemperatureLevelOffset t1(pressure_sigma, temp_expect, 0, true);
  sv.add_observer(t1);

  for(int i = 0; i < temp_expect.rows(); ++i) {
    BOOST_CHECK_CLOSE(t1.temperature(pressure_sigma->pressure_grid()(i)).convert(units::K).value.value(), temp_expect(i), 1e-3);
  }
}

BOOST_AUTO_TEST_SUITE_END()


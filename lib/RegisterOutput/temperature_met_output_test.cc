#include "unit_test_support.h"
#include "met_data_fixture.h"
#include "output_hdf.h"
#include "temperature_met_output.h"

using namespace FullPhysics;
using namespace blitz;
BOOST_FIXTURE_TEST_SUITE(temperature_met_output, MetDataFixture)

BOOST_AUTO_TEST_CASE(basic)
{
  TemperatureMetOutput po(boost::shared_ptr<TemperatureMet>(new TemperatureMet(met_data, pressure, 0, true)));
  boost::shared_ptr<OutputHdf> out(new OutputHdf("temperature_met_output.h5", 20, 112, 5, 3));
  add_file_to_cleanup("temperature_met_output.h5");
  po.register_output(out);

  // Simple test, we just make sure that we can write output. All the
  // actual value calculation is checked in temperature unit test.

  out->write();
  BOOST_CHECK(true);
}

BOOST_AUTO_TEST_SUITE_END()

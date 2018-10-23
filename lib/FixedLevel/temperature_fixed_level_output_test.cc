#include "unit_test_support.h"
#include "temperature_fixed_level_output.h"
#include "output_hdf.h"
#include "met_data_fixture.h"

using namespace FullPhysics;
using namespace blitz;
BOOST_FIXTURE_TEST_SUITE(temperature_fixed_level_output, MetDataFixture)

BOOST_AUTO_TEST_CASE(basic)
{
    boost::shared_ptr<PressureLevelInput> press_level(new PressureLevelInput(met_data->pressure_levels()));
    Array<double, 1> temp = met_data->temperature();
    Array<bool, 1> flag_temp(temp.rows());
    flag_temp = false;
    bool flag_offset = false;
    double toffset = 0;
    boost::shared_ptr<TemperatureFixedLevel> temp_met(new TemperatureFixedLevel(flag_temp, flag_offset, temp, toffset, pressure, press_level));

    TemperatureFixedLevelOutput po(temp_met);
    boost::shared_ptr<OutputHdf> out(new OutputHdf("temperature_fixed_level_output.h5", 20, 112, 5, 3));
    add_file_to_cleanup("temperature_fixed_level_output.h5");
    po.register_output(out);

    // Simple test, we just make sure that we can write output. All the
    // actual value calculation is checked in temperature unit test.

    out->write();
    BOOST_CHECK(true);
}

BOOST_AUTO_TEST_SUITE_END()



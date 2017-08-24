#include "temperature_fixed_level.h"
#include "met_data_fixture.h"
#include "unit_test_support.h"
#include "pressure_fixed_level.h"

using namespace FullPhysics;
using namespace blitz;

BOOST_FIXTURE_TEST_SUITE(temperature_fixed_level, MetDataFixture)

BOOST_AUTO_TEST_CASE(basic)
{
    boost::shared_ptr<PressureLevelInput> press_level(new PressureLevelInput(met_data->pressure_levels()));
    StateVector sv;
    Array<bool, 1> flag_temp(3);
    Array<double, 1> temp(3);
    bool flag_offset;
    flag_temp= true, true, false;
    flag_offset = false;
    temp = 10, 11, 12;
    double toffset = 0;
    TemperatureFixedLevel t1(flag_temp, flag_offset, temp, toffset, pressure, press_level);
    sv.add_observer(t1);
    Array<double, 1> temp_expect(3);
    temp_expect = 10, 11, 12;
    BOOST_CHECK_MATRIX_CLOSE(t1.temperature_levels().value(), temp_expect);
    Array<double, 1> x(2);
    x = 1, 2;
    sv.update_state(x);
    temp_expect = 1, 2, 12;
    BOOST_CHECK_MATRIX_CLOSE(t1.temperature_levels().value(), temp_expect);
    
    flag_temp = false, false, false;
    flag_offset = 1;
    StateVector sv2;
    temp = 1, 2, 12;
    TemperatureFixedLevel t2(flag_temp, flag_offset, temp, toffset, pressure, press_level);
    sv2.add_observer(t2);
    temp_expect = 1, 2, 12;
    BOOST_CHECK_MATRIX_CLOSE(t2.temperature_levels().value(), temp_expect);
    sv2.update_state(x);
    temp_expect = 1 + 1, 2 + 1, 12 + 1;
    BOOST_CHECK_MATRIX_CLOSE(t2.temperature_levels().value(), temp_expect);
}

BOOST_AUTO_TEST_CASE(met)
{
    boost::shared_ptr<PressureLevelInput> press_level(new PressureLevelInput(met_data->pressure_levels()));
    Array<double, 1> temp = met_data->temperature();
    Array<bool, 1> flag_temp(temp.rows());
    flag_temp = false;
    bool flag_offset = false;
    double toffset = 0;
    TemperatureFixedLevel temp_met(flag_temp, flag_offset, temp, toffset, pressure, press_level);

    IfstreamCs expected_data(test_data_dir() + "expected/absorber/ecmwf");
    Array<double, 1> texpect(19);
    texpect =
        205.821929931640625, 215.999053955078125, 226.8797454833984375, 241.315185546875,
        254.6123046875, 262.12969970703125, 262.55548095703125, 256.668182373046875,
        246.0504608154296875, 236.890838623046875, 230.515411376953125, 224.625244140625,
        220.57806396484375, 216.279083251953125, 213.21875, 210.87371826171875,
        210.322174072265625, 210.211090087890625, 160.11948565552248169;

    Array<double, 1> tcalc = temp_met.temperature_grid(*pressure).value.value();
    BOOST_CHECK_MATRIX_CLOSE(tcalc, texpect);
}

BOOST_AUTO_TEST_SUITE_END()

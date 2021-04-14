#include "temperature_level.h"
#include "unit_test_support.h"
#include "met_data_fixture.h"

using namespace FullPhysics;
using namespace blitz;

// Consolidate reused test setup and checking
class TemperatureLevelTest {
public:
    TemperatureLevelTest(const boost::shared_ptr<Pressure>& press_in) : pressure(press_in) {
        test_temp.resize(19);
        test_temp = 244.2, 214.553, 218.029, 222.544, 218.341, 221.37, 227.38,
            233.493, 239.376, 244.52, 248.708, 251.979, 254.537, 256.655, 258.521,
            260.155, 261.747, 261.732, 258.598;

        temp_lev.reset(new TemperatureLevel(test_temp, pressure));
    }

    void check_results(const boost::shared_ptr<TemperatureLevel>& check_obj) {
        for(int i = 0; i < test_temp.rows() - 1; ++i) {
            // Check interpolation by simply picking a pressure between two levels for
            // a simple check
            AutoDerivativeWithUnit<double> grid_pressure_1 = pressure->pressure_grid()(i);
            AutoDerivativeWithUnit<double> grid_pressure_2 = pressure->pressure_grid()(i+1);
            AutoDerivativeWithUnit<double> check_pressure((grid_pressure_1.value + grid_pressure_2.value)/2, grid_pressure_1.units);
            double expt_temperature = (test_temp(i) + test_temp(i+1)) / 2;

            BOOST_CHECK_CLOSE(check_obj->temperature(check_pressure).convert(units::K).value.value(), expt_temperature, 1e-8);
        }
    }

    Array<double, 1> test_temp;
    boost::shared_ptr<Pressure> pressure;
    boost::shared_ptr<TemperatureLevel> temp_lev;
};

BOOST_FIXTURE_TEST_SUITE(temperature_level, MetDataFixture)

BOOST_AUTO_TEST_CASE(basic)
{

    TemperatureLevelTest temp_test(pressure);

    StateVector sv;
    sv.add_observer(*temp_test.temp_lev);

    temp_test.check_results(temp_test.temp_lev);
}

BOOST_AUTO_TEST_CASE(serialization)
{
    TemperatureLevelTest temp_test(pressure);

    std::string serial_str = serialize_write_string(temp_test.temp_lev);
    boost::shared_ptr<TemperatureLevel> temp_lev_read = serialize_read_string<TemperatureLevel>(serial_str);

    temp_test.check_results(temp_lev_read);
}

BOOST_AUTO_TEST_SUITE_END()


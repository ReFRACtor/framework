#include "unit_test_support.h"
#include "hdf_file.h"
#include "example_met_file.h"

using namespace FullPhysics;
using namespace blitz;

BOOST_FIXTURE_TEST_SUITE(example_met_file, GlobalFixture)

BOOST_AUTO_TEST_CASE(load_file)
{
    boost::shared_ptr<HdfFile> h_file(new HdfFile(test_data_dir() + "in/meteorology/example_met_data.h5"));
    ExampleMetFile met_file(h_file, "20091009203401");

    blitz::Array<double, 1> pressure_levels_read = h_file->read_field<double, 2>("/Meteorology/pressure_levels")(0, Range::all());
    BOOST_CHECK_MATRIX_CLOSE_TOL(met_file.pressure_levels(), pressure_levels_read, 1e-12);

    blitz::Array<double, 1> specific_humidity_read = h_file->read_field<double, 2>("/Meteorology/specific_humidity")(0, Range::all());
    BOOST_CHECK_MATRIX_CLOSE_TOL(met_file.specific_humidity(), specific_humidity_read, 1e-12);

    blitz::Array<double, 1> temperature_read = h_file->read_field<double, 2>("/Meteorology/temperature")(0, Range::all());
    BOOST_CHECK_MATRIX_CLOSE_TOL(met_file.temperature(), temperature_read, 1e-12);

    double surface_pressure_read = h_file->read_field<double, 1>("/Meteorology/surface_pressure")(0);
    BOOST_CHECK_CLOSE(met_file.surface_pressure(), surface_pressure_read, 1e-12);

    double windspeed_u_read = h_file->read_field<double, 1>("/Meteorology/windspeed_u")(0);
    BOOST_CHECK_CLOSE(met_file.windspeed_u(), windspeed_u_read, 1e-12);

    double windspeed_v_read = h_file->read_field<double, 1>("/Meteorology/windspeed_v")(0);
    BOOST_CHECK_CLOSE(met_file.windspeed_v(), windspeed_v_read, 1e-12);

}

BOOST_AUTO_TEST_SUITE_END()

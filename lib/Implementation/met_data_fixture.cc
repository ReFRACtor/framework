#include "met_data_fixture.h"
#include "hdf_file.h"
#include "pressure_sigma.h"
#include "example_met_file.h"

using namespace FullPhysics;
using namespace blitz;

MetDataFixture::MetDataFixture()
{
    IfstreamCs press_data(test_data_dir() + "/in/meteorology/pressure_in");

    blitz::Array<double, 1> pressure_array;
    double psurf_val;

    press_data >> pressure_array;
    psurf_val = pressure_array(pressure_array.rows()-1);

    boost::shared_ptr<HdfFile> met_file(new HdfFile(test_data_dir() + "in/meteorology/example_met_data.h5"));
    met_data.reset(new ExampleMetFile(met_file, "20091009203401"));

    pressure.reset(new PressureSigma(pressure_array, psurf_val));
}

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
}

BOOST_AUTO_TEST_SUITE_END()

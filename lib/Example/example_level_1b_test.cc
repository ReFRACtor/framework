#include "unit_test_support.h"
#include "hdf_file.h"
#include "example_level_1b.h"

using namespace FullPhysics;
using namespace blitz;

BOOST_FIXTURE_TEST_SUITE(example_level1b, GlobalFixture)

BOOST_AUTO_TEST_CASE(load_file)
{
    boost::shared_ptr<HdfFile> h_file(new HdfFile(test_data_dir() + "expected/level_1b_example/example_l1b.h5"));
    ExampleLevel1b l1b_file(h_file, "20091009203401");
}

BOOST_AUTO_TEST_SUITE_END()

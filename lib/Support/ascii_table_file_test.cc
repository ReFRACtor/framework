#include "ascii_table_file.h"
#include "global_fixture.h"
#include "unit_test_support.h"

using namespace FullPhysics;
using namespace blitz;

BOOST_FIXTURE_TEST_SUITE(ascii_table_file, GlobalFixture)

BOOST_AUTO_TEST_CASE(simple)
{
    // Test that we can read a simple ascii file
    AsciiTableFile ascii_file(input_dir() + "cross_sections/no2r_97.nm");

    Array<double, 2> data(ascii_file.data());

    // Check first and last row as well as sizes
    BOOST_CHECK_EQUAL(data.rows(), 21485);
    BOOST_CHECK_EQUAL(data.cols(), 2);

    BOOST_CHECK_CLOSE(data(0, 0), 279.9995873, 1e-8);
    BOOST_CHECK_CLOSE(data(0, 1), 5.40253E-20, 1e-8);

    BOOST_CHECK_CLOSE(data(data.rows()-1, 0), 666.7621421, 1e-8);
    BOOST_CHECK_CLOSE(data(data.rows()-1, 1), 8.82147E-21, 1e-8);
}

BOOST_AUTO_TEST_CASE(comments)
{
    // Test that we can read a simple ascii file
    AsciiTableFile ascii_file(input_dir() + "cross_sections/o3abs_brion_195_660_vacfinal.dat");

    Array<double, 2> data(ascii_file.data());

    // Check first and last row as well as sizes
    BOOST_CHECK_EQUAL(data.rows(), 46501);
    BOOST_CHECK_EQUAL(data.cols(), 4);

    BOOST_CHECK_CLOSE(data(0, 0), 195.0640, 1e-8);
    BOOST_CHECK_CLOSE(data(0, 1), 3.8231E+01, 1e-8);
    BOOST_CHECK_CLOSE(data(0, 2), 8.7657E-03, 1e-8);
    BOOST_CHECK_CLOSE(data(0, 3), 1.5841E-04, 1e-8);

    // Last row with non-zero columns 2 and 3
    BOOST_CHECK_CLOSE(data(32400, 0), 519.1445, 1e-8);
    BOOST_CHECK_CLOSE(data(32400, 1), 1.7713E-01, 1e-8);
    BOOST_CHECK_CLOSE(data(32400, 2), 5.5255E-05, 1e-8);
    BOOST_CHECK_CLOSE(data(32400, 3), 7.9594E-07, 1e-8);

    BOOST_CHECK_CLOSE(data(data.rows()-1, 0), 660.1823, 1e-8);
    BOOST_CHECK_CLOSE(data(data.rows()-1, 1), 2.1001E-01, 1e-8);
    BOOST_CHECK_CLOSE(data(data.rows()-1, 2), 0, 1e-8);
    BOOST_CHECK_CLOSE(data(data.rows()-1, 3), 0, 1e-8);
}

BOOST_AUTO_TEST_SUITE_END()

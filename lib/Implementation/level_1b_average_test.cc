#include "unit_test_support.h"
#include "level_1b_average.h"
#include "example_level_1b.h"
#include "fp_exception.h"

using namespace FullPhysics;
using namespace blitz;

BOOST_FIXTURE_TEST_SUITE(level_1b_average, GlobalFixture)

BOOST_AUTO_TEST_CASE(basic)
{
    boost::shared_ptr<HdfFile> hf(new HdfFile(test_data_dir() + "in/common/l1b_example_data.h5"));

    std::vector<boost::shared_ptr<Level1b> > l1b;

    boost::shared_ptr<Level1b> l1b_1(new ExampleLevel1b(hf, "2014090915251774"));
    l1b.push_back(l1b_1);

    boost::shared_ptr<Level1b> l1b_2(new ExampleLevel1b(hf, "2014120112331638"));
    l1b.push_back(l1b_2);

    Level1bAverage h(l1b);
    BOOST_CHECK_EQUAL(h.number_spectrometer(), 3);
    blitz::Array<double, 2> st_expect(3,4);
    st_expect = 
        0.5, 0.3677656278, -0.2218676396, 0.0,
        0.5, 0.3678720668, -0.2218151002, 0.0,
        0.5, 0.3678224012, -0.2217771686, 0.0;

    for(int i = 0; i < 3; ++i) {
        // Check that these are close enough, the value per spectrometer does vary
        BOOST_CHECK_CLOSE(h.latitude(i).value, -3.850e+01, 1e-1);
        BOOST_CHECK_CLOSE(h.solar_zenith(i).value, 3.986e+01, 1e-1);
        BOOST_CHECK_CLOSE(h.solar_azimuth(i).value, 3.1112e+02, 1e-2);
        BOOST_CHECK_CLOSE(h.altitude(i).value, 0.0, 1e-2);
        BOOST_CHECK_CLOSE(h.sounding_zenith(i).value, 3.148e+01, 1e-1);
        BOOST_CHECK_CLOSE(h.sounding_azimuth(i).value, 1.303e+02, 1e-1);

        BOOST_CHECK_MATRIX_CLOSE_TOL(h.stokes_coefficient(i), st_expect(i, Range::all()), 1e-4);
    }
    BOOST_CHECK_CLOSE(h.relative_velocity(0).value, 3.3503947754e+03, 1e-4);
    BOOST_CHECK_CLOSE(h.time(0).pgs_time(), 6.8442992565e+08, 1e-4);

    IfstreamCs expected_data(test_data_dir() + "expected/level_1b_average/basic");
    Array<double, 1> rad_expected, uncert_expected;
    expected_data >> rad_expected;
    expected_data >> uncert_expected;
    BOOST_CHECK_MATRIX_CLOSE_TOL(h.radiance(0).data(), rad_expected, 1e10);
    BOOST_CHECK_MATRIX_CLOSE_TOL(h.radiance(0).uncertainty(), uncert_expected, 1e10);
}

BOOST_AUTO_TEST_SUITE_END()

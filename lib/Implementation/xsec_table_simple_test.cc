#include "xsec_table_simple.h"
#include "global_fixture.h"
#include "unit_test_support.h"
#include "fp_serialize_support.h"
#include "ascii_table_file.h"

using namespace FullPhysics;
using namespace blitz;

BOOST_FIXTURE_TEST_SUITE(xsec_table_simple, GlobalFixture)

const boost::shared_ptr<XSecTableSimple> construct_xsec_simple(const std::string& filename) {
    // Open cross section table file and read into appropriate values
    AsciiTableFile xsec_table_file(filename);

    Array<double, 2> xsec_table_values(xsec_table_file.data());

    ArrayWithUnit<double, 1> spectral_grid(xsec_table_values(Range::all(), 0), Unit("nm"));
    ArrayWithUnit<double, 2> xsec_values(xsec_table_values(Range::all(), Range(1,1)), Unit("cm^2"));

    double conversion_factor = 1;

    boost::shared_ptr<XSecTableSimple> xsec_simple(new XSecTableSimple(spectral_grid, xsec_values, conversion_factor));

    return xsec_simple;
}

void check_xsec_simple(const boost::shared_ptr<XSecTableSimple>& xsec_simple, const std::string& expt_filename)
{
    AsciiTableFile expected_interp(expt_filename);
    Array<double, 2> expected_data(expected_interp.data());

    double sum_sq = 0;
    for(int grid_idx = 0; grid_idx < expected_data.rows(); grid_idx++) {
        DoubleWithUnit spectral_point(expected_data(grid_idx, 0), Unit("nm"));
        DoubleWithUnit xsec_value = xsec_simple->cross_section_value(spectral_point);
        
        // Scale up to get a more comparable RMS value
        sum_sq += pow(expected_data(grid_idx, 1)*1e20 - xsec_value.value*1e20, 2);

        BOOST_CHECK_CLOSE(expected_data(grid_idx, 1), xsec_value.value, 1.4e-1);
    }

    double rms = sqrt(sum_sq / expected_data.rows());

    // Check that RMS difference mathches this threshold
    BOOST_CHECK_LT(rms, 4e-3);
}

BOOST_AUTO_TEST_CASE(basic)
{
    boost::shared_ptr<XSecTableSimple> xsec_simple = construct_xsec_simple(input_dir() + "cross_sections/no2r_97.nm");
    check_xsec_simple(xsec_simple, test_data_dir() + "expected/cross_section/no2_xsec_interp_value.dat");
}

BOOST_AUTO_TEST_CASE(serialization)
{
    boost::shared_ptr<XSecTableSimple> xsec_orig = construct_xsec_simple(input_dir() + "cross_sections/no2r_97.nm");
    std::string serial_str = serialize_write_string(xsec_orig);
    boost::shared_ptr<XSecTableSimple> xsec_reconstructed = serialize_read_string<XSecTableSimple>(serial_str);
    check_xsec_simple(xsec_reconstructed, test_data_dir() + "expected/cross_section/no2_xsec_interp_value.dat");
}

BOOST_AUTO_TEST_SUITE_END()

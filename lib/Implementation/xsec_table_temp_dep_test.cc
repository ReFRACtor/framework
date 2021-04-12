#include "xsec_table_temp_dep.h"
#include "global_fixture.h"
#include "unit_test_support.h"
#include "fp_serialize_support.h"
#include "ascii_table_file.h"

using namespace FullPhysics;
using namespace blitz;

BOOST_FIXTURE_TEST_SUITE(xsec_table_temp_dep, GlobalFixture)

const boost::shared_ptr<XSecTableTempDep> construct_xsec_temp_dep(const std::string& filename) {
    // Open cross section table file and read into appropriate values
    AsciiTableFile xsec_table_file(filename);

    Array<double, 2> xsec_table_values(xsec_table_file.data());

    ArrayWithUnit<double, 1> spectral_grid(xsec_table_values(Range::all(), 0), Unit("nm"));
    ArrayWithUnit<double, 2> xsec_values(xsec_table_values(Range::all(), Range(1,3)), Unit("cm^2"));

    double conversion_factor = 1e20;

    boost::shared_ptr<XSecTableTempDep> xsec_temp_dep(new XSecTableTempDep(spectral_grid, xsec_values, conversion_factor));

    return xsec_temp_dep;
}

void check_xsec_temp_dep(const boost::shared_ptr<XSecTableTempDep>& xsec_temp_dep, const std::string& expt_filename)
{
    AsciiTableFile expected_interp(expt_filename);
    Array<double, 2> expected_data(expected_interp.data());

    double sum_sq = 0;
    for(int grid_idx = 0; grid_idx < expected_data.rows(); grid_idx++) {
        DoubleWithUnit spectral_point(expected_data(grid_idx, 0), Unit("nm"));
        DoubleWithUnit xsec_value = xsec_temp_dep->cross_section_value(spectral_point);
        
        // Scale up to get a more comparable RMS value
        sum_sq += pow(expected_data(grid_idx, 1)*1e20 - xsec_value.value*1e20, 2);

        BOOST_CHECK_CLOSE(expected_data(grid_idx, 1), xsec_value.value, 8.8e-2);
    }

    double rms = sqrt(sum_sq / expected_data.rows());

    // Check that RMS difference mathches this threshold
    BOOST_CHECK_LT(rms, 8e-3);
}

BOOST_AUTO_TEST_CASE(basic)
{
    boost::shared_ptr<XSecTableTempDep> xsec_temp_dep = construct_xsec_temp_dep(input_dir() + "cross_sections/o3abs_brion_195_660_vacfinal.dat");
    check_xsec_temp_dep(xsec_temp_dep, test_data_dir() + "expected/cross_section/o3_xsec_interp_value.dat");
}

BOOST_AUTO_TEST_CASE(serialization)
{
    boost::shared_ptr<XSecTableTempDep> xsec_orig = construct_xsec_temp_dep(input_dir() + "cross_sections/o3abs_brion_195_660_vacfinal.dat");
    std::string serial_str = serialize_write_string(xsec_orig);
    boost::shared_ptr<XSecTableTempDep> xsec_reconstructed = serialize_read_string<XSecTableTempDep>(serial_str);
    check_xsec_temp_dep(xsec_reconstructed, test_data_dir() + "expected/cross_section/o3_xsec_interp_value.dat");
}

BOOST_AUTO_TEST_SUITE_END()

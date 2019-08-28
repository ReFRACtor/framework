#include "unit_test_support.h"
#include "global_fixture.h"

#include "hdf_file.h"

#include "pca_optical_properties.h"
#include "pca_binning.h"

using namespace FullPhysics;
using namespace blitz;

BOOST_FIXTURE_TEST_SUITE(pca_binning, GlobalFixture)

BOOST_AUTO_TEST_CASE(compare_with_offline)
{
    // Contains optical data inputs as well as expected outputs in the format
    // produced in the offline Fortran code
    // Does not contain all the values accessible by PCAOpticalProperties interface
    boost::shared_ptr<HdfFile> test_data(new HdfFile(test_data_dir() + "in/pca/binning_test_data.h5"));
    boost::shared_ptr<PCAOpticalProperties> opt_props(new PCAOpticalPropertiesFile(test_data));

    // Matches 
    auto num_bin_points_shape = test_data->read_shape<1>("num_bin_points");
    int num_bins = num_bin_points_shape(0);
    auto binning = PCABinning(opt_props, num_bins);

    // Get expected values
    auto num_bin_points_expt = test_data->read_field<int, 1>("num_bin_points");
    auto bin_indexes_expt_packed = test_data->read_field<int, 1>("bin_indexes_packed");

    // Compare number of bin points, if this isn't the same the following certainly won't work
    std::cerr << "expt = " << num_bin_points_expt << std::endl;
    std::cerr << "calc = " << binning.num_bin_points() << std::endl;
    BOOST_CHECK_MATRIX_CLOSE(num_bin_points_expt, binning.num_bin_points());

    // Compare bin indexes while unpacking the saved values
    int pack_start = 0;
    for(int bidx = 0; bidx < num_bins; bidx++) {
        int npoints = num_bin_points_expt(bidx);

        // Remove 1 from expected values to make them zero based indexes
        BOOST_CHECK_MATRIX_CLOSE(bin_indexes_expt_packed(Range(pack_start, pack_start + npoints-1)) - 1, binning.bin_indexes()[bidx]);
        pack_start += npoints;
    }
}

BOOST_AUTO_TEST_SUITE_END()

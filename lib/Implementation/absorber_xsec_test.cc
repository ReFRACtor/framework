#include "unit_test_support.h"
#include "global_fixture.h"
#include "fp_serialize_support.h"

#include "absorber_xsec.h"
#include "atmosphere_standard.h"
#include "ascii_table_file.h"

using namespace FullPhysics;
using namespace blitz;

BOOST_FIXTURE_TEST_SUITE(absorber_xsec, GlobalFixture)

BOOST_AUTO_TEST_CASE(basic)
{
    // This unit test relies upon serialized data created by the Python creator system
    // The serialized data is recreated by this script:
    // test/unit/serialization/serialize_uv_atmosphere.py
    // Changes that do not impact the stored data do not require serialization to be rerun.
    // But if data variables are added or the type of dimensions change then the serialization will be rerun.
    // Remember that the data is the only thing serialized, the current code is always rerun. Serialization
    // just makes it easier to set up a class hiearchy with canned data.
    //
    // The expected data here comes from the uvvis_od_calc repository
#ifdef FP_HAVE_BOOST_SERIALIZATION
    std::string atmosphere_serialized_fn = test_data_dir() + "in/uv_atmosphere/uv_atmosphere.xml";
    boost::shared_ptr<AtmosphereStandard> uv_atmosphere = serialize_read<AtmosphereStandard>(atmosphere_serialized_fn);
    boost::shared_ptr<AbsorberXSec> abs_xsec = boost::dynamic_pointer_cast<AbsorberXSec>(uv_atmosphere->absorber_ptr());

    // Preload the expected data, the first column is the grid, rest of columns are value per layer
    std::vector<Array<double, 2> > expected_data;

    // Expected values from from the uvvis_od_calc repository
    int expt_number_cols = -1;
    for(int gas_idx = 0; gas_idx < abs_xsec->number_species(); gas_idx++) {
        std::string gas_name = abs_xsec->gas_name(gas_idx);
        std::string gas_expt_fn = test_data_dir() + "expected/cross_section/absorber_od_" + gas_name + ".dat";
        AsciiTableFile expt_file(gas_expt_fn);
        expected_data.push_back(expt_file.data());

        // Check that all data has same number of columns
        if(expt_number_cols < 0) {
            expt_number_cols = expt_file.data().cols();
        } else if (expt_number_cols != expt_file.data().cols()) {
            Exception err;
            err << "File " << gas_expt_fn << " has " << expt_file.data().cols() << " columns instead of the expected " << expt_number_cols;
            throw err;
        }
    }

    // Each expected file contains the grid, pull the grid from the first file and first column
    Array<double, 1> grid_nm(expected_data[0](Range::all(), 0));
    Range expt_data_cols(1, expt_number_cols-1);

    if(grid_nm.rows() == 0 || expt_number_cols == 0) {
        throw Exception("Error loading expected results data");
    }

    // Only 1 sensor defined in uv_atmosphere_config.py
    int sensor_idx = 0;

    //for(int grid_idx = 0; grid_idx < grid_nm.rows(); grid_idx++) {
    { int grid_idx = 0;
        float wn = DoubleWithUnit(grid_nm(grid_idx), Unit("nm")).convert_wave(units::inv_cm).value;
        ArrayAdWithUnit<double, 2> comp_abs_od = abs_xsec->optical_depth_each_layer(wn, sensor_idx);

        for(int gas_idx = 0; gas_idx < abs_xsec->number_species(); gas_idx++) {
            Array<double, 1> gas_expt_abs_od(expected_data[gas_idx](grid_idx, expt_data_cols));
            Array<double, 1> gas_comp_abs_od(comp_abs_od.value.value()(Range::all(), gas_idx));
            BOOST_CHECK_MATRIX_CLOSE_TOL(gas_expt_abs_od, gas_comp_abs_od, 1.5e-5);
        }
    }
#else
    throw Exception("Can not run this unit test without serialization support");
#endif
}

BOOST_AUTO_TEST_SUITE_END()


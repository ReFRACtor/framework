#include "unit_test_support.h"
#include "global_fixture.h"
#include "fp_serialize_support.h"

#include "absorber_xsec.h"
#include "atmosphere_standard.h"
#include "ascii_table_file.h"

using namespace FullPhysics;
using namespace blitz;

BOOST_FIXTURE_TEST_SUITE(absorber_xsec, GlobalFixture)

BOOST_AUTO_TEST_CASE(optical_depth)
{
  // This unit test relies upon serialized data created by the Python
  // creator system The serialized data is recreated by this script:
  // test/unit/serialization/serialize_uv_atmosphere.py Changes that
  // do not impact the stored data do not require serialization to be
  // rerun.  But if data variables are added or the type of dimensions
  // change then the serialization will be rerun.  Remember that the
  // data is the only thing serialized, the current code is always
  // rerun. Serialization just makes it easier to set up a class
  // hiearchy with canned data.
  //
  // The expected data here comes from the uvvis_od_calc repository
  if(!have_serialize_supported())
    return;
  
  std::string atmosphere_serialized_fn = test_data_dir() + "in/uv_atmosphere/uv_atmosphere.bin.gz";
  boost::shared_ptr<AtmosphereStandard> uv_atmosphere = serialize_read_binary<AtmosphereStandard>(atmosphere_serialized_fn);
  boost::shared_ptr<AbsorberXSec> abs_xsec = boost::dynamic_pointer_cast<AbsorberXSec>(uv_atmosphere->absorber_ptr());
  
  // Only 1 sensor defined in uv_atmosphere_config.py
  int sensor_idx = 0;
  
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
  
  for(int grid_idx = 0; grid_idx < grid_nm.rows(); grid_idx++) {
    float wn = DoubleWithUnit(grid_nm(grid_idx), Unit("nm")).convert_wave(units::inv_cm).value;
    ArrayAdWithUnit<double, 2> comp_abs_od = abs_xsec->optical_depth_each_layer(wn, sensor_idx);
    
    for(int gas_idx = 0; gas_idx < abs_xsec->number_species(); gas_idx++) {
      Array<double, 1> gas_expt_abs_od(expected_data[gas_idx](grid_idx, expt_data_cols));
      Array<double, 1> gas_comp_abs_od(comp_abs_od.value.value()(Range::all(), gas_idx));
      BOOST_CHECK_MATRIX_CLOSE_TOL(gas_expt_abs_od, gas_comp_abs_od, 6e-5);
    }
  }
}

BOOST_AUTO_TEST_CASE(density)
{
  if(!have_serialize_supported())
    return;

  std::string atmosphere_serialized_fn = test_data_dir() + "in/uv_atmosphere/uv_atmosphere.bin.gz";
  boost::shared_ptr<AtmosphereStandard> uv_atmosphere = serialize_read_binary<AtmosphereStandard>(atmosphere_serialized_fn);
  boost::shared_ptr<AbsorberXSec> abs_xsec = boost::dynamic_pointer_cast<AbsorberXSec>(uv_atmosphere->absorber_ptr());
  
  // Compare total air density against offline code, make sure we have commensurate units by converting our results
  // to those of the expected values
  std::string total_density_expt_fn = test_data_dir() + "expected/cross_section/total_air_density.dat";
  AsciiTableFile total_air_expt_file(total_density_expt_fn);
  
  // Comparison data in the second column
  Array<double, 1> total_air_expt(total_air_expt_file.data()(Range::all(), 1));
  Array<double, 1> total_air_calc(abs_xsec->total_air_number_density_level().convert(Unit("cm^-3")).value.value());
  
  // Normalize values so we are not comparing large numbers
  total_air_expt = total_air_expt / max(total_air_expt);
  total_air_calc = total_air_calc / max(total_air_calc);
  
  BOOST_CHECK_MATRIX_CLOSE_TOL(total_air_expt, total_air_calc, 1e-8);
  
  // Check per gas density values
  Array<double, 2> gas_density_expt(abs_xsec->number_layer() + 1, abs_xsec->number_species());
  for(int gas_idx = 0; gas_idx < abs_xsec->number_species(); gas_idx++) {
    std::string gas_name = abs_xsec->gas_name(gas_idx);
    std::string gas_expt_fn = test_data_dir() + "expected/cross_section/gas_air_density_" + gas_name + ".dat";
    AsciiTableFile expt_file(gas_expt_fn);
    
    // Expected values in second column
    gas_density_expt(Range::all(), gas_idx) = expt_file.data()(Range::all(), 1);
  }
  
  Array<double, 2> gas_density_calc(abs_xsec->gas_number_density_level().convert(Unit("cm^-3")).value.value());
  
  // Normalize values so we are not comparing large numbers
  gas_density_expt = gas_density_expt / max(gas_density_expt);
  gas_density_calc = gas_density_calc / max(gas_density_calc);
  
  BOOST_CHECK_MATRIX_CLOSE_TOL(gas_density_expt, gas_density_calc, 1e-8);
}



BOOST_AUTO_TEST_SUITE_END()


#include "absorber_absco.h"
#include "lua_configuration_fixture.h"
#include "unit_test_support.h"
#include "atmosphere_fixture.h"
#include "lua_state.h"
#include <cstdlib>

using namespace FullPhysics;
using namespace blitz;

BOOST_FIXTURE_TEST_SUITE(absorber_absco, LuaConfigurationFixture)

BOOST_AUTO_TEST_CASE(optical_depth_each_layer)
{
  secondIndex i2;
  IfstreamCs expected_data(test_data_dir() + "expected/absorber_absco/optical_depth_each_layer");
  Absorber& a = *config_absorber;
  Array<double, 1> od_expect;
  expected_data >> od_expect;
  Array<double, 1> od_calc_1( sum(a.optical_depth_each_layer(12929.94,0).value(),i2) );
  BOOST_CHECK_MATRIX_CLOSE_TOL(od_calc_1, od_expect, 1e-10);
  expected_data >> od_expect;
  Array<double, 1> od_calc_2( sum(a.optical_depth_each_layer(12930.30,0).value(),i2) );
  BOOST_CHECK_MATRIX_CLOSE_TOL(od_calc_2, od_expect, 1e-10);
  if (false) {
    std::cerr << setprecision(20) << std::scientific
              << "# od_expect 12929.94" << std::endl
              << od_calc_1 << std::endl
              << "# od_expect 12930.30" << std::endl
              << od_calc_2 << std::endl;
  }
}

BOOST_AUTO_TEST_CASE(optical_depth_each_layer_direct_integrate)
{
  secondIndex i2;
  IfstreamCs expected_data(test_data_dir() + 
                   "expected/absorber_absco/optical_depth_each_layer");
  AbsorberAbsco& a = dynamic_cast<AbsorberAbsco&>(*config_absorber);
  Array<double, 1> od_expect;
  expected_data >> od_expect;
  Array<double, 1> od_calc_1( sum(a.optical_depth_each_layer_direct_integrate(12929.94,0),i2) );
  BOOST_CHECK_MATRIX_CLOSE_TOL(od_calc_1, od_expect, 1e-8);
  expected_data >> od_expect;
  Array<double, 1> od_calc_2( sum(a.optical_depth_each_layer_direct_integrate(12930.30,0),i2) );
  BOOST_CHECK_MATRIX_CLOSE_TOL(od_calc_2, od_expect, 1e-8);
}

BOOST_AUTO_TEST_CASE(press_wf)
{
  IfstreamCs expected_data(test_data_dir() + 
                           "expected/absorber_absco/pressure_wf");
  AbsorberAbsco& a = dynamic_cast<AbsorberAbsco&>(*config_absorber);
  Array<double, 1> press_wf_lay_expect, press_wf_lev_expect;
  expected_data >> press_wf_lay_expect >> press_wf_lev_expect;
  BOOST_CHECK_MATRIX_CLOSE(a.pressure_weighting_function_layer(0).value(), 
                           press_wf_lay_expect);
  BOOST_CHECK_MATRIX_CLOSE(a.pressure_weighting_function_grid(0).value(), 
                           press_wf_lev_expect);
  BOOST_CHECK_CLOSE(a.xgas(0, "CO2").value(), 0.00039670981755915669, 1e-4);
  if (false) {
    std::cerr << setprecision(20) << std::scientific
              << "# Expected pressure weighting function for layers." << std::endl
              << a.pressure_weighting_function_layer(0).value() << std::endl
              << "# Expected pressure weighting function for levels" << std::endl
              << a.pressure_weighting_function_grid(0).value() << std::endl;
  }    
}

BOOST_AUTO_TEST_CASE(gas_column_thickness)
{
  // Make sure that wet - dry column == h2o column as it should
  AbsorberAbsco& a = dynamic_cast<AbsorberAbsco&>(*config_absorber);
  Array<double, 1> wet_col(a.wet_air_number_density_layer(0).value.value());
  Array<double, 1> dry_col(a.dry_air_number_density_layer(0).value.value());
  Array<double, 1> h2o_col(a.gas_column_thickness_layer(0, "H2O").value.value());
  // The two comparison values are on the order of 1e23 so must have a higher tol here
  BOOST_CHECK_MATRIX_CLOSE_TOL(wet_col - dry_col, h2o_col, 1e13);

  // Check that dry / column(O2) = vmr
  Array<double, 1> o2_col(a.gas_column_thickness_layer(0, "O2").value.value());
  Array<double, 1> o2_vmr(o2_col.rows());
  o2_vmr = a.average_vmr("O2").value();
  BOOST_CHECK_MATRIX_CLOSE(o2_col / dry_col, o2_vmr);
}

BOOST_AUTO_TEST_CASE(jacobian)
{
  AbsorberAbsco& a = dynamic_cast<AbsorberAbsco&>(*config_absorber);
  StateVector& sv = *config_state_vector;
  Array<double, 1> sv0(sv.state().copy());
  // Pick an band where we have some CO2, so varying the VMR of CO2 affects 
  // the results.
  int spec_index = 2;
  double wn = 4820.0;
  ArrayAd<double, 2> od = a.optical_depth_each_layer(wn, spec_index);
  Array<double, 2> od0(od.shape());
  od0 = od.value();
  Array<double, 3> jac = od.jacobian().copy();
  for(int i = 0; i < sv.state().rows(); ++i) {
    Array<double, 1> svn(sv0.copy());
    svn(i) += epsilon(i);
    sv.update_state(svn);
    Array<double, 2> jacfd(od0.shape());
    jacfd = (a.optical_depth_each_layer_nder(wn,spec_index) - od0) 
      / epsilon(i);
    if(false) {                        // Can turn this off to dump values,
                                // if needed for debugging
      double diff = max(abs(jac(Range::all(), Range::all(), i) - jacfd));
      if(diff > 0)
        std::cerr << i << ": " << jac(Range::all(), Range::all(), i) << "\n"
                  << jacfd << "\n"
                  << diff << "\n";
    }
    BOOST_CHECK_MATRIX_CLOSE_TOL(jac(Range::all(), Range::all(), i), jacfd, 1e-8);
  }
}

BOOST_AUTO_TEST_CASE(serialization)
{
  secondIndex i2;
  if(!have_serialize_supported())
    return;
  
  AbsorberAbsco& a = dynamic_cast<AbsorberAbsco&>(*config_absorber);
  serialize_write("absorber_absco.xml", config_absorber);
  std::string d = serialize_write_string(config_absorber);
  if(false)
    std::cerr << d;
  boost::shared_ptr<AbsorberAbsco> ar =
    serialize_read_string<AbsorberAbsco>(d);
  Array<double, 1> od_calc_1_1( sum(a.optical_depth_each_layer(12929.94,0).
				    value(),i2) );
  Array<double, 1> od_calc_1_2( sum(ar->optical_depth_each_layer(12929.94,0).
				    value(),i2) );
  Array<double, 1> od_calc_1_3( sum(a.optical_depth_each_layer_direct_integrate(12929.94,0),i2) );
  Array<double, 1> od_calc_1_4( sum(ar->optical_depth_each_layer_direct_integrate(12929.94,0),i2) );
  BOOST_CHECK_MATRIX_CLOSE_TOL(od_calc_1_1, od_calc_1_2, 1e-10);
  BOOST_CHECK_MATRIX_CLOSE_TOL(od_calc_1_1, od_calc_1_3, 1e-10);
  BOOST_CHECK_MATRIX_CLOSE_TOL(od_calc_1_1, od_calc_1_4, 1e-10);
  Array<double, 1> od_calc_2_1( sum(a.optical_depth_each_layer(12930.30,0).
				  value(),i2) );
  Array<double, 1> od_calc_2_2( sum(a.optical_depth_each_layer(12930.30,0).
				  value(),i2) );
  Array<double, 1> od_calc_2_3( sum(a.optical_depth_each_layer_direct_integrate(12930.30,0),i2) );
  Array<double, 1> od_calc_2_4( sum(ar->optical_depth_each_layer_direct_integrate(12930.30,0),i2) );
  BOOST_CHECK_MATRIX_CLOSE_TOL(od_calc_2_1, od_calc_2_2, 1e-10);
  BOOST_CHECK_MATRIX_CLOSE_TOL(od_calc_2_1, od_calc_2_3, 1e-10);
  BOOST_CHECK_MATRIX_CLOSE_TOL(od_calc_2_1, od_calc_2_4, 1e-10);
  BOOST_CHECK_MATRIX_CLOSE(a.pressure_weighting_function_layer(0).value(), 
			   ar->pressure_weighting_function_layer(0).value());
  BOOST_CHECK_MATRIX_CLOSE(a.pressure_weighting_function_grid(0).value(), 
			   ar->pressure_weighting_function_grid(0).value());

  BOOST_CHECK_MATRIX_CLOSE(a.wet_air_number_density_layer(0).value.value(),
			   ar->wet_air_number_density_layer(0).value.value());
  BOOST_CHECK_MATRIX_CLOSE(a.dry_air_number_density_layer(0).value.value(),
			   ar->dry_air_number_density_layer(0).value.value());
   BOOST_CHECK_MATRIX_CLOSE(a.total_air_number_density_layer(0).value.value(),
			   ar->total_air_number_density_layer(0).value.value());
  BOOST_CHECK_MATRIX_CLOSE(a.gas_column_thickness_layer(0, "H2O").value.value(),
			   ar->gas_column_thickness_layer(0, "H2O").value.value());
  BOOST_CHECK_MATRIX_CLOSE(a.gas_column_thickness_layer(0, "O2").value.value(),
			   ar->gas_column_thickness_layer(0, "O2").value.value());
  // Pick an band where we have some CO2, so varying the VMR of CO2 affects 
  // the results.
  int spec_index = 2;
  double wn = 4820.0;
  ArrayAd<double, 2> od_1 = a.optical_depth_each_layer(wn, spec_index);
  ArrayAd<double, 2> od_2 = ar->optical_depth_each_layer(wn, spec_index);
  BOOST_CHECK_MATRIX_CLOSE(od_1.jacobian(), od_2.jacobian());
}
BOOST_AUTO_TEST_SUITE_END()

BOOST_FIXTURE_TEST_SUITE(absorber_absco2, ConfigurationTwoBroadener)
BOOST_AUTO_TEST_CASE(two_broadener_optical_depth_each_layer)
{
  secondIndex i2;
  // This include a two broadener case, a 1 broadener case, and no
  // broadener case to check the logic for all of these.

  // Test 2 broadener
  AbsorberAbsco& a = *(boost::dynamic_pointer_cast<AbsorberAbsco>(config_absorber));
  Array<double, 1> od_calc_1( sum(a.optical_depth_each_layer(6200.10,0).value(),i2) );
  Array<double, 1> od_calc_2( sum(a.optical_depth_each_layer_nder(6200.10,0),i2) );
  std::cerr << od_calc_1;
  std::cerr << od_calc_2;
  // Test one broadener
  Array<double, 1> od_calc_3( sum(a.optical_depth_each_layer(4820.0,2).value(),i2) );
  Array<double, 1> od_calc_4( sum(a.optical_depth_each_layer_nder(4820.0,2),i2) );
  std::cerr << od_calc_3;
  std::cerr << od_calc_4;
}
BOOST_AUTO_TEST_SUITE_END()

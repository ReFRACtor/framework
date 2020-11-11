#include "ils_table.h"
#include "unit_test_support.h"

using namespace FullPhysics;
using namespace blitz;

BOOST_FIXTURE_TEST_SUITE(ils_table, GlobalFixture)

BOOST_AUTO_TEST_CASE(table)
{
  HdfFile hf(test_data_dir() + "in/ils/ils_linear_table.h5");
  IlsTableLinear ils(hf, 0, "A-Band", "o2");
  BOOST_CHECK_EQUAL(ils.band_name(), "A-Band");
  BOOST_CHECK_EQUAL(ils.hdf_band_name(), "o2");
  IfstreamCs expected(test_data_dir() + "expected/ils_table/table");
  double wn_center_val;
  Array<double, 1> wn, response_expected;
  expected >> wn_center_val >> wn >> response_expected;
  Array<double, 1> grad(2);
  grad = 1,0.1;
  AutoDerivative<double> wn_center(wn_center_val, grad);
  ArrayAd<double, 1> response;
  ils.ils(wn_center, wn, response);
  BOOST_CHECK_MATRIX_CLOSE(response.value(), response_expected);
  double epsilon = 0.001;
  ArrayAd<double, 1> response2;
  ils.ils(wn_center_val + epsilon, wn, response2);
  BOOST_CHECK_MATRIX_CLOSE(response.jacobian()(Range::all(), 0),
                           (response2.value() - response.value()) / epsilon);
  BOOST_CHECK_MATRIX_CLOSE(response.jacobian()(Range::all(), 1),
                           (response2.value() - response.value()) / epsilon * 0.1);
}

// Repeat tests from above, but with the function_type set to interpol
BOOST_AUTO_TEST_CASE(interpol)
{
  HdfFile hf(test_data_dir() + "in/ils/ils_interpol.h5");
  IlsTableLinear ils(hf, 0, "A-Band", "o2");
  BOOST_CHECK_EQUAL(ils.band_name(), "A-Band");
  BOOST_CHECK_EQUAL(ils.hdf_band_name(), "o2");
  IfstreamCs expected(test_data_dir() + "expected/ils_table/interpol");
  double wn_center_val;
  Array<double, 1> wn, response_expected;
  expected >> wn_center_val >> wn >> response_expected;
  Array<double, 1> grad(2);
  grad = 1,0.1;
  AutoDerivative<double> wn_center(wn_center_val, grad);
  ArrayAd<double, 1> response;
  ils.ils(wn_center, wn, response);
  BOOST_CHECK_MATRIX_CLOSE(response.value(), response_expected);
  double epsilon = 0.001;
  ArrayAd<double, 1> response2;
  ils.ils(wn_center_val + epsilon, wn, response2);
  // Typical value of jacobian is 0.001, so 1e-6 is about 0.1% difference
  // with finite difference which is perfectly reasonable.
  BOOST_CHECK_MATRIX_CLOSE_TOL(response.jacobian()(Range::all(), 0),
          (response2.value() - response.value()) / epsilon, 1e-6);
  BOOST_CHECK_MATRIX_CLOSE_TOL(response.jacobian()(Range::all(), 1),
          (response2.value() - response.value()) / epsilon * 0.1, 1e-6);
}

// Make sure can use a table where values are interpreted in log space
BOOST_AUTO_TEST_CASE(log_table)
{
  // Use OCO ILS values
  HdfFile hf(test_data_dir() + "in/ils/ils_log_table.h5");
  IlsTableLog ils(hf, 0, "A-Band", "o2");
  BOOST_CHECK_EQUAL(ils.band_name(), "A-Band");
  BOOST_CHECK_EQUAL(ils.hdf_band_name(), "o2");
  IfstreamCs expected(test_data_dir() + "expected/ils_table/log_table");
  double center_val;
  Array<double, 1> wl, response_expected;
  expected >> center_val >> wl >> response_expected;
  Array<double, 1> grad(2);
  grad = 1,0.1;
  AutoDerivative<double> wl_center(center_val, grad);
  ArrayAd<double, 1> response;
  ils.ils(wl_center, wl, response);
  BOOST_CHECK_MATRIX_CLOSE_TOL(response.value(), response_expected, 1e-4);
}

BOOST_AUTO_TEST_CASE(serialization)
{
  if(!have_serialize_supported())
    return;
  HdfFile hf(test_data_dir() + "in/ils/ils_linear_table.h5");
  boost::shared_ptr<IlsTableLinear> ils =
    boost::make_shared<IlsTableLinear>(hf, 0, "A-Band", "o2");
  std::string d = serialize_write_string(ils);
  if(false)
    std::cerr << d;
  boost::shared_ptr<IlsTableLinear> ils_r = serialize_read_string<IlsTableLinear>(d);
  BOOST_CHECK_EQUAL(ils_r->band_name(), "A-Band");
  BOOST_CHECK_EQUAL(ils_r->hdf_band_name(), "o2");
  IfstreamCs expected(test_data_dir() + "expected/ils_table/table");
  double wn_center_val;
  Array<double, 1> wn, response_expected;
  expected >> wn_center_val >> wn >> response_expected;
  Array<double, 1> grad(2);
  grad = 1,0.1;
  AutoDerivative<double> wn_center(wn_center_val, grad);
  ArrayAd<double, 1> response;
  ils_r->ils(wn_center, wn, response);
  BOOST_CHECK_MATRIX_CLOSE(response.value(), response_expected);
  double epsilon = 0.001;
  ArrayAd<double, 1> response2;
  ils_r->ils(wn_center_val + epsilon, wn, response2);
  BOOST_CHECK_MATRIX_CLOSE(response.jacobian()(Range::all(), 0),
                           (response2.value() - response.value()) / epsilon);
  BOOST_CHECK_MATRIX_CLOSE(response.jacobian()(Range::all(), 1),
                           (response2.value() - response.value()) / epsilon * 0.1);
}
BOOST_AUTO_TEST_SUITE_END()

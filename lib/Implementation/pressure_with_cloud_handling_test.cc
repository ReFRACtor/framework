#include "pressure_with_cloud_handling.h"
#include "pressure_sigma.h"
#include "unit_test_support.h"
#include "fp_serialize_support.h"

using namespace FullPhysics;
using namespace blitz;

class PrintNotify: public Observer<Pressure> {
public:
  PrintNotify() {}
  virtual ~PrintNotify() {}
  virtual void notify_update(const Pressure& P)
  {
    std::cerr << "Notify_update: " << P << "\n";
  }
};
  
BOOST_FIXTURE_TEST_SUITE(pressure_with_cloud_handling, GlobalFixture)

BOOST_AUTO_TEST_CASE(basic)
{
  Array<double, 1> a(3), b(3);
  a = 0; b = 0.3, 0.6, 1.0;
  double psurf = 10;
  double cloud_pressure_level = 7;
  auto psigma = boost::make_shared<PressureSigma>(a, b, psurf);
  auto psigma2 = boost::make_shared<PressureSigma>(a.reverse(blitz::firstDim),
                                                   b.reverse(blitz::firstDim),
                                                   psurf, Pressure::PREFER_DECREASING_PRESSURE);
  PressureWithCloudHandling p(psigma, cloud_pressure_level);
  PressureWithCloudHandling p2(psigma2, cloud_pressure_level);
  BOOST_CHECK_EQUAL(p.type_preference(), Pressure::PREFER_INCREASING_PRESSURE);
  BOOST_CHECK_EQUAL(p2.type_preference(), Pressure::PREFER_DECREASING_PRESSURE);
  PrintNotify obs;
  p.add_observer(obs);
  Array<double, 1> press_grid_expect(3);
  press_grid_expect = 3, 6, 10;
  Array<double, 1> press_grid_expect2(3);
  press_grid_expect2 = 10, 6, 3;
  BOOST_CHECK_MATRIX_CLOSE(p.pressure_grid().value.value(), press_grid_expect);
  BOOST_CHECK_MATRIX_CLOSE(p.pressure_grid(Pressure::NATIVE_ORDER).value.value(), press_grid_expect);
  BOOST_CHECK_MATRIX_CLOSE(p2.pressure_grid().value.value(), press_grid_expect);
  BOOST_CHECK_MATRIX_CLOSE(p2.pressure_grid(Pressure::NATIVE_ORDER).value.value(), press_grid_expect2);
  p.do_cloud(true);
  p2.do_cloud(true);
  Array<double, 1> press_grid_expect3(2);
  press_grid_expect3 = 3, 6;
  Array<double, 1> press_grid_expect4(2);
  press_grid_expect4 = 6, 3;
  BOOST_CHECK_MATRIX_CLOSE(p.pressure_grid().value.value(), press_grid_expect3);
  BOOST_CHECK_MATRIX_CLOSE(p.pressure_grid(Pressure::NATIVE_ORDER).value.value(), press_grid_expect3);
  BOOST_CHECK_MATRIX_CLOSE(p2.pressure_grid().value.value(), press_grid_expect3);
  BOOST_CHECK_MATRIX_CLOSE(p2.pressure_grid(Pressure::NATIVE_ORDER).value.value(), press_grid_expect4);
  psigma->set_surface_pressure(15);
  cloud_pressure_level = 0.0001;
  PressureWithCloudHandling p3(psigma, cloud_pressure_level);
  p3.do_cloud(true);
  BOOST_CHECK_THROW(p3.pressure_grid(), Exception);
  BOOST_CHECK_THROW(p3.pressure_grid(Pressure::DECREASING_PRESSURE), Exception);
}

BOOST_AUTO_TEST_CASE(serialization)
{
  if(!have_serialize_supported())
    return;
  Array<double, 1> a(3), b(3);
  a = 0; b = 0.3, 0.6, 1.0;
  double psurf = 10;
  double cloud_pressure_level = 7;
  auto p = boost::make_shared<PressureWithCloudHandling>
    (boost::make_shared<PressureSigma>(a,b, psurf), cloud_pressure_level);
  std::string d = serialize_write_string(p);
  auto pr = serialize_read_string<PressureWithCloudHandling>(d);
  Array<double, 1> press_grid_expect(3);
  press_grid_expect = 3, 6, 10;
  BOOST_CHECK_MATRIX_CLOSE(pr->pressure_grid().value.value(), press_grid_expect);
  pr->do_cloud(true);
  Array<double, 1> press_grid_expect2(2);
  press_grid_expect2 = 3, 6;
  BOOST_CHECK_MATRIX_CLOSE(pr->pressure_grid().value.value(), press_grid_expect2);
}

BOOST_AUTO_TEST_SUITE_END()


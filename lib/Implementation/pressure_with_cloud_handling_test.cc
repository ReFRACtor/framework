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
  auto psigma = boost::make_shared<PressureSigma>(a, b, psurf);
  PressureWithCloudHandling p(psigma, 7);
  PrintNotify obs;
  p.add_observer(obs);
  Array<double, 1> press_grid_expect(3);
  press_grid_expect = 3, 6, 10;
  BOOST_CHECK_MATRIX_CLOSE(p.pressure_grid().value.value(), press_grid_expect);
  std::cerr << p.pressure_grid().value.value() << "\n";
  p.do_cloud(true);
  std::cerr << p.pressure_grid().value.value() << "\n";
  Array<double, 1> press_grid_expect2(2);
  press_grid_expect2 = 3, 6;
  BOOST_CHECK_MATRIX_CLOSE(p.pressure_grid().value.value(), press_grid_expect2);
  psigma->set_surface_pressure(15);
}

BOOST_AUTO_TEST_CASE(serialization)
{
  if(!have_serialize_supported())
    return;
  Array<double, 1> a(3), b(3);
  a = 0; b = 0.3, 0.6, 1.0;
  double psurf = 10;
  auto p = boost::make_shared<PressureWithCloudHandling>
    (boost::make_shared<PressureSigma>(a,b, psurf), 7);
  std::string d = serialize_write_string(p);
  if(false)
    std::cerr << d;
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

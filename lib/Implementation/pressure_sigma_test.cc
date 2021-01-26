#include "pressure_sigma.h"
#include "configuration_fixture.h"
#include "unit_test_support.h"
#include "fp_serialize_support.h"
#include "generic_object_map.h"

using namespace FullPhysics;
using namespace blitz;

// Test class to make sure that Observer are properly serialized.
namespace FullPhysics {
class TestObserver : public Observer<Pressure> {
public:
  TestObserver() : data(0) {}
  virtual ~TestObserver() {}
  virtual void notify_update(const Pressure& P)
  { data = P.surface_pressure_value(); };
  int data;
  template<class Archive>
  void serialize(Archive & ar, const unsigned int)
  {
    ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(ObserverPressure)
      & FP_NVP(data);
  }
};
}

FP_EXPORT_KEY(TestObserver);

FP_IMPLEMENT(TestObserver);

BOOST_FIXTURE_TEST_SUITE(pressure_sigma, GlobalFixture)

BOOST_AUTO_TEST_CASE(basic)
{
  StateVector sv;
  Array<double, 1> a(3), b(3);
  a = 0; b = 0.3, 0.6, 1.0;
  double psurf = 10;
  PressureSigma p(a,b, psurf);
  Array<double, 1> press_grid_expect(3);
  press_grid_expect = 3, 6, 10;
  BOOST_CHECK_CLOSE(p.surface_pressure().value.value(), psurf, 1e-4);
  BOOST_CHECK_MATRIX_CLOSE(p.pressure_grid().value.value(), press_grid_expect);
  sv.add_observer(p);
  Array<double, 1> x(1);
  x = 20;
  sv.update_state(x);
  BOOST_CHECK_CLOSE(p.surface_pressure().value.value(), 20, 1e-4);
  press_grid_expect = 6, 12, 20;
  BOOST_CHECK_MATRIX_CLOSE(p.pressure_grid().value.value(), press_grid_expect);
}

BOOST_AUTO_TEST_CASE(serialization)
{
  if(!have_serialize_supported())
    return;
  Array<double, 1> a(3), b(3);
  a = 0; b = 0.3, 0.6, 1.0;
  double psurf = 10;
  boost::shared_ptr<PressureSigma> p = boost::make_shared<PressureSigma>(a,b, psurf);
  boost::shared_ptr<TestObserver> pobs = boost::make_shared<TestObserver>();
  p->add_observer(*pobs);
  std::string d = serialize_write_string(p);
  if(false)
    std::cerr << d;
  boost::shared_ptr<PressureSigma> pr = serialize_read_string<PressureSigma>(d);
  Array<double, 1> press_grid_expect(3);
  press_grid_expect = 3, 6, 10;
  BOOST_CHECK_CLOSE(pr->surface_pressure().value.value(), psurf, 1e-4);
  BOOST_CHECK_MATRIX_CLOSE(pr->pressure_grid().value.value(), press_grid_expect);
  // Since TestObserver wasn't serialized, shouldn't have a connection
  // with pr.
  pr->set_surface_pressure(100);
  BOOST_CHECK(fabs(pobs->data - pr->surface_pressure_value()) > 1e-4);
}

BOOST_AUTO_TEST_CASE(serialization2)
{
  // Second serialization, where we have an observer that is serialized.
  if(!have_serialize_supported())
    return;
  Array<double, 1> a(3), b(3);
  a = 0; b = 0.3, 0.6, 1.0;
  double psurf = 10;
  boost::shared_ptr<GenericObjectMap> m =
    boost::make_shared<GenericObjectMap>();
  (*m)["p"] = boost::make_shared<PressureSigma>(a, b, psurf);
  (*m)["pobs"] = boost::make_shared<TestObserver>();
  m->get<PressureSigma>("p")->add_observer(*m->get<TestObserver>("pobs"));
  std::string d = serialize_write_string(m);
  if(false)
    std::cerr << d;
  boost::shared_ptr<GenericObjectMap> mr =
    serialize_read_string<GenericObjectMap>(d);
  Array<double, 1> press_grid_expect(3);
  press_grid_expect = 3, 6, 10;
  BOOST_CHECK_CLOSE(mr->get<PressureSigma>("p")->surface_pressure().value.value(), psurf, 1e-4);
  BOOST_CHECK_MATRIX_CLOSE(mr->get<PressureSigma>("p")->pressure_grid().value.value(), press_grid_expect);
  // Since TestObserver was serialized, we should have a connection
  // with tsr->p.
  mr->get<PressureSigma>("p")->set_surface_pressure(100);
  BOOST_CHECK(fabs(mr->get<TestObserver>("pobs")->data -
	   mr->get<PressureSigma>("p")->surface_pressure_value()) < 1e-4);
}

// BOOST_AUTO_TEST_CASE(jacobian)
// {
//   Pressure& p = *config_pressure;
//   StateVector& sv = *config_state_vector;
//   ArrayAd<double, 1> pgrid = p.pressure_grid();
//   Array<double, 1> pgrid0(pgrid.value().copy());
//   Array<double, 2> jac = pgrid.jacobian();
//   Array<double, 1> sv0(sv.state().copy());
//   for(int i = 0; i < sv.state().rows(); ++i) {
//     Array<double, 1> svn(sv0.copy());
//     svn(i) += epsilon(i);
//     sv.update_state(svn);
//     Array<double, 1> jacfd(pgrid0.shape());
//     jacfd = (p.pressure_grid().value() - pgrid0) / epsilon(i);
//     if(false) {			// Can turn this off to dump values,
// 				// if needed for debugging
//       double diff = max(abs(jac(Range::all(), i) - jacfd));
//       if(diff > 0)
// 	std::cerr << i << ": " << jac(Range::all(), i) << "\n"
// 		  << jacfd << "\n"
// 		  << diff << "\n";
//     }
//     BOOST_CHECK_MATRIX_CLOSE_TOL(jac(Range::all(), i), jacfd, 1e-3);
//   }
// }

BOOST_AUTO_TEST_SUITE_END()


#include "unit_test_support.h"
#include "l_rad_rt.h"
#include "lsi_rt.h"
#include "lidort_fixture.h"

using namespace FullPhysics;
using namespace blitz;

// Perhaps move this elsewhere
class LRadLambertianFixture : public LidortLambertianFixture {
public:
  LRadLambertianFixture()
  {
    // Set viewing geometry
    sza = 74.128288268999995;
    zen = 30.0;
    azm = 10.0;

    rt_first_order.reset(new LRadRt(lidort_rt, config_spectral_window->spectral_bound(), 
                                    sza, zen, azm, pure_nadir, true, false));
    rt_second_order.reset(new LRadRt(lidort_rt, config_spectral_window->spectral_bound(),
                                     sza, zen, azm, pure_nadir, true, true));
    rt_lrad_only.reset(new LRadRt(lidort_rt->stokes_coefficient(),
                                  config_atmosphere,
                                  config_spectral_window->spectral_bound(),
                                  sza, zen, azm, pure_nadir,
                                  lidort_rt->number_stokes()));
    rt_lrad_only_second.reset(new LRadRt(lidort_rt->stokes_coefficient(),
                                  config_atmosphere,
                                  config_spectral_window->spectral_bound(),
                                  sza, zen, azm, pure_nadir,
                                  lidort_rt->number_stokes(), true));
    check_tol = 1e-6;
  }
  virtual ~LRadLambertianFixture() {}

  double check_tol;
  boost::shared_ptr<LRadRt> rt_first_order;
  boost::shared_ptr<LRadRt> rt_second_order;
  boost::shared_ptr<LRadRt> rt_lrad_only;
  boost::shared_ptr<LRadRt> rt_lrad_only_second;
};

class LRadCoxmunkFixture : public LidortCoxmunkFixture {
public:
  LRadCoxmunkFixture()
  {
    // Set viewing geometry
    sza = 74.128288268999995;
    zen = 30.0;
    azm = 10.0;

    rt_first_order.reset(new LRadRt(lidort_rt, config_spectral_window->spectral_bound(), 
                                        sza, zen, azm, pure_nadir, true, false));
    rt_second_order.reset(new LRadRt(lidort_rt, config_spectral_window->spectral_bound(),
                                         sza, zen, azm, pure_nadir, true, true));
    rt_lrad_only.reset(new LRadRt(lidort_rt->stokes_coefficient(),
                                      config_atmosphere,
                                      config_spectral_window->spectral_bound(),
                                      sza, zen, azm, pure_nadir,
                                      lidort_rt->number_stokes()));
    rt_lrad_only_second.reset(new LRadRt(lidort_rt->stokes_coefficient(),
                                      config_atmosphere,
                                      config_spectral_window->spectral_bound(),
                                      sza, zen, azm, pure_nadir,
                                      lidort_rt->number_stokes(), true));
    check_tol = 1e-6;
  }
  virtual ~LRadCoxmunkFixture() {}

  double check_tol;
  boost::shared_ptr<LRadRt> rt_first_order;
  boost::shared_ptr<LRadRt> rt_second_order;
  boost::shared_ptr<LRadRt> rt_lrad_only;
  boost::shared_ptr<LRadRt> rt_lrad_only_second;
};

// --------- lambertian ----------- //

BOOST_FIXTURE_TEST_SUITE(l_rad_rt_lambertian, LRadLambertianFixture)

BOOST_AUTO_TEST_CASE(reflectance_first_order)
{
  Array<double, 1> wn(2);
  wn = 12929.94, 12930.30;
  Array<double, 2> stokes_expect(2, 3);
  stokes_expect = 
    0.020742234580895982377, -0.0013228620534916765358, -0.00047622341054522510529, 
    0.020742000112342170309, -0.001322989829824955171, -0.0004762694093559555219;
  BOOST_CHECK_MATRIX_CLOSE_TOL(rt_first_order->stokes(wn, 0),
                               stokes_expect, check_tol);
  BOOST_CHECK_MATRIX_CLOSE_TOL
    (rt_first_order->stokes_and_jacobian(wn, 0).value(),
     stokes_expect, check_tol);
  Array<double, 1> rad_expect(stokes_expect(Range::all(), 0));
  BOOST_CHECK_MATRIX_CLOSE_TOL(rt_first_order->reflectance(wn, 0, true).
                               spectral_range().data(), 
                               rad_expect, check_tol);
  BOOST_CHECK_MATRIX_CLOSE_TOL
    (rt_first_order->reflectance(wn, 0).spectral_range().data(), 
     rad_expect, check_tol);
}

BOOST_AUTO_TEST_CASE(reflectance_first_order_serialization)
{
  if(!have_serialize_supported())
    return;
  std::string d = serialize_write_string(rt_first_order);
  if(false)
    std::cerr << d;
  boost::shared_ptr<LRadRt> rt_first_order_r =
    serialize_read_string<LRadRt>(d);
  Array<double, 1> wn(2);
  wn = 12929.94, 12930.30;
  Array<double, 2> stokes_expect(2, 3);
  stokes_expect = 
    0.020742234580895982377, -0.0013228620534916765358, -0.00047622341054522510529, 
    0.020742000112342170309, -0.001322989829824955171, -0.0004762694093559555219;
  BOOST_CHECK_MATRIX_CLOSE_TOL(rt_first_order_r->stokes(wn, 0),
                               stokes_expect, check_tol);
  BOOST_CHECK_MATRIX_CLOSE_TOL
    (rt_first_order_r->stokes_and_jacobian(wn, 0).value(),
     stokes_expect, check_tol);
  Array<double, 1> rad_expect(stokes_expect(Range::all(), 0));
  BOOST_CHECK_MATRIX_CLOSE_TOL(rt_first_order_r->reflectance(wn, 0, true).
                               spectral_range().data(), 
                               rad_expect, check_tol);
  BOOST_CHECK_MATRIX_CLOSE_TOL
    (rt_first_order_r->reflectance(wn, 0).spectral_range().data(), 
     rad_expect, check_tol);
}

BOOST_AUTO_TEST_CASE(reflectance_lrad_only)
{
  Array<double, 1> wn(2);
  wn = 12929.94, 12930.30;
  Array<double, 2> stokes_expect(2, 3);
  stokes_expect = 
    0.0094660708716226022591, -0.0013228620534916765358, -0.00047622341054522510529, 
    0.0094659533466465214241, -0.001322989829824955171, -0.0004762694093559555219;
  BOOST_CHECK_MATRIX_CLOSE_TOL(rt_lrad_only->stokes(wn, 0),
                               stokes_expect, check_tol);
  BOOST_CHECK_MATRIX_CLOSE_TOL(rt_lrad_only->stokes_and_jacobian(wn, 0).value(),
                               stokes_expect, check_tol);
}

BOOST_AUTO_TEST_CASE(reflectance_second_order)
{
  Array<double, 1> wn(2);
  wn = 12929.94, 12930.30;
  Array<double, 2> stokes_expect(2, 3);
  stokes_expect = 
    0.020725928098217986928, -0.001440501649682172057, -0.00052790819320807464907, 
    0.02072569062801354195, -0.0014406412567393635932, -0.00052795921740695468064;
  BOOST_CHECK_MATRIX_CLOSE_TOL(rt_second_order->stokes(wn, 0), 
                               stokes_expect, check_tol);
  BOOST_CHECK_MATRIX_CLOSE_TOL
    (rt_second_order->stokes_and_jacobian(wn, 0).value(), 
     stokes_expect, check_tol);
  Array<double, 1> rad_expect(stokes_expect(Range::all(), 0));
  BOOST_CHECK_MATRIX_CLOSE_TOL(rt_second_order->reflectance(wn, 0, true).
                               spectral_range().data(), 
                               rad_expect, check_tol);
  BOOST_CHECK_MATRIX_CLOSE_TOL
    (rt_second_order->reflectance(wn, 0).spectral_range().data(), 
     rad_expect, check_tol);
}

BOOST_AUTO_TEST_CASE(reflectance_lrad_only_serialization)
{
  if(!have_serialize_supported())
    return;
  std::string d = serialize_write_string(rt_lrad_only);
  if(false)
    std::cerr << d;
  boost::shared_ptr<LRadRt> rt_lrad_only_r =
    serialize_read_string<LRadRt>(d);
  Array<double, 1> wn(2);
  wn = 12929.94, 12930.30;
  Array<double, 2> stokes_expect(2, 3);
  stokes_expect = 
    0.0094660708716226022591, -0.0013228620534916765358, -0.00047622341054522510529, 
    0.0094659533466465214241, -0.001322989829824955171, -0.0004762694093559555219;
  BOOST_CHECK_MATRIX_CLOSE_TOL(rt_lrad_only_r->stokes(wn, 0),
                               stokes_expect, check_tol);
  BOOST_CHECK_MATRIX_CLOSE_TOL(rt_lrad_only_r->stokes_and_jacobian(wn, 0).value(),
                               stokes_expect, check_tol);
}

BOOST_AUTO_TEST_CASE(reflectance_second_order_serialization)
{
  if(!have_serialize_supported())
    return;
  std::string d = serialize_write_string(rt_second_order);
  if(false)
    std::cerr << d;
  boost::shared_ptr<LRadRt> rt_second_order_r =
    serialize_read_string<LRadRt>(d);
  Array<double, 1> wn(2);
  wn = 12929.94, 12930.30;
  Array<double, 2> stokes_expect(2, 3);
  stokes_expect = 
    0.020725928098217986928, -0.001440501649682172057, -0.00052790819320807464907, 
    0.02072569062801354195, -0.0014406412567393635932, -0.00052795921740695468064;
  BOOST_CHECK_MATRIX_CLOSE_TOL(rt_second_order_r->stokes(wn, 0), 
                               stokes_expect, check_tol);
  BOOST_CHECK_MATRIX_CLOSE_TOL
    (rt_second_order_r->stokes_and_jacobian(wn, 0).value(), 
     stokes_expect, check_tol);
  Array<double, 1> rad_expect(stokes_expect(Range::all(), 0));
  BOOST_CHECK_MATRIX_CLOSE_TOL(rt_second_order_r->reflectance(wn, 0, true).
                               spectral_range().data(), 
                               rad_expect, check_tol);
  BOOST_CHECK_MATRIX_CLOSE_TOL
    (rt_second_order_r->reflectance(wn, 0).spectral_range().data(), 
     rad_expect, check_tol);
}

BOOST_AUTO_TEST_CASE(l_rad_timing)
{
  is_timing_test();
  boost::timer tm;
  RtAtmosphere& atm = *config_atmosphere;
  int i = 0;
  for(double wn = 12929.94; wn <= 13210.15; wn += 0.01) {
    ArrayAd<double, 1> 
      res(rt_lrad_only->stokes_and_jacobian_single_wn(wn, 0));
    if(++i % 1000 == 0)
      std::cerr << "Done with " << i << "\n"
                << "Total: " << tm.elapsed() << "\n"
                << atm.timer_info() << "\n";
  }
}

BOOST_AUTO_TEST_CASE(jac_first_order)
{
  // We use rt_lrad_only so we only test the portion of the Jacobian coming
  // from LRad. We separately test LIDORT.
  is_long_test();
  StateVector& sv = *config_state_vector;
  Array<double, 1> sv0(config_initial_guess->initial_guess());
  sv.update_state(sv0);
  int spec_index = 2;
  Array<double, 1> wn(1);
  wn = 4820.0;

  // Test jacobians for all three l_rad_first pseudo spherical modes
  for(int ps_idx = LRadDriver::REGULAR; ps_idx <= LRadDriver::PLANE_PARALLEL; ps_idx++) {
    LRadDriver::PsMode ps_mode = static_cast<LRadDriver::PsMode>(ps_idx);
    boost::shared_ptr<LRadRt> lrad_ps(new LRadRt(rt_lrad_only->stokes_coefficient(),
                                                         config_atmosphere,
                                                         config_spectral_window->spectral_bound(),
                                                         sza, zen, azm, pure_nadir,
                                                         rt_lrad_only->number_stokes(),
                                                         false, 4, 0.01,
                                                         ps_mode));
    ArrayAd<double, 2> stk = 
      lrad_ps->stokes_and_jacobian(wn, spec_index);
    Array<double, 2> stk0(stk.shape());
    stk0 = stk.value();
    Array<double, 3> jac = stk.jacobian().copy();
    for(int i = 0; i < sv.state().rows(); ++i) {
      Array<double, 1> svn(sv0.copy());
      svn(i) += epsilon(i);
      sv.update_state(svn);
      Array<double, 2> jacfd(stk.shape());
      jacfd = (lrad_ps->stokes(wn, spec_index) - stk0) / epsilon(i);
      double diff = max(abs(jac(Range::all(), Range::all(), i) - jacfd));
      if(false) {                        // Can turn this on to dump values,
                                  // if needed for debugging
        if(diff > 0) {
          std::cerr << i << ": " << diff << "\n"
                    << diff / max(abs(jacfd)) << "\n"
                    << jac(0,Range::all(), i) << "\n"
                    << jacfd << "\n";
        }
      }
      if(diff > 1e-6) {
        std::cerr << i << ": " << diff << "\n"
                  << diff / max(abs(jacfd)) << "\n"
                  << jac(0,Range::all(), i) << "\n"
                  << jacfd << "\n";
      }
      BOOST_CHECK(diff < 1e-6);
    }
  }
}

BOOST_AUTO_TEST_CASE(jac_second_order)
{
  // We use rt_lrad_only so we only test the portion of the Jacobian coming
  // from LRad. We separately test LIDORT.
  is_long_test();
  StateVector& sv = *config_state_vector;
  Array<double, 1> sv0(config_initial_guess->initial_guess());
  sv.update_state(sv0);
  int spec_index = 2;
  Array<double, 1> wn(1);
  wn = 4820.0;
  ArrayAd<double, 2> stk = 
    rt_lrad_only_second->stokes_and_jacobian(wn, spec_index);
  Array<double, 2> stk0(stk.shape());
  stk0 = stk.value();
  Array<double, 3> jac = stk.jacobian().copy();
  for(int i = 0; i < sv.state().rows(); ++i) {
    Array<double, 1> svn(sv0.copy());
    svn(i) += epsilon(i);
    sv.update_state(svn);
    Array<double, 2> jacfd(stk.shape());
    jacfd = (rt_lrad_only_second->stokes(wn, spec_index) - stk0) 
      / epsilon(i);
    double diff = max(abs(jac(Range::all(), Range::all(), i) - jacfd));
    if(false) {                        // Can turn this on to dump values,
                                // if needed for debugging
      if(diff > 0) {
        std::cerr << i << ": " << diff << "\n"
                  << diff / max(abs(jacfd)) << "\n"
                  << jac(0,Range::all(), i) << "\n"
                  << jacfd << "\n";
      }
    }
    if(diff > 1e-4) {
      std::cerr << i << ": " << diff << "\n"
                << diff / max(abs(jacfd)) << "\n"
                << jac(0,Range::all(), i) << "\n"
                << jacfd << "\n";
    }
    BOOST_CHECK(diff < 1e-4);
  }
}

BOOST_AUTO_TEST_CASE(regular_ps)
{
  // Set up to run l_rad in regular ps mode
  Array<double, 1> zen_small(3);
  zen_small = 1.0e-6;
  boost::shared_ptr<LRadRt> lrad_ps(new LRadRt(rt_lrad_only->stokes_coefficient(),
                                                        config_atmosphere,
                                                        config_spectral_window->spectral_bound(),
                                                        sza, zen_small, azm, pure_nadir,
                                                        rt_lrad_only->number_stokes(),
                                                        false, 4, 0.01,
                                                        LRadDriver::REGULAR));
  Array<double, 1> wn(1);
  wn = 12929.94;
  Array<double, 1> stokes_expect(3);
  stokes_expect = 
    0.0088645271123479803949, -0.0011307484227398416094, -0.00041155876643655497481;
  BOOST_CHECK_MATRIX_CLOSE_TOL(lrad_ps->stokes(wn, 0)(0, Range::all()), stokes_expect, 1e-6);
}

BOOST_AUTO_TEST_CASE(regular_ps_serialization)
{
  // Set up to run l_rad in regular ps mode
  if(!have_serialize_supported())
    return;
  Array<double, 1> zen_small(3);
  zen_small = 1.0e-6;
  boost::shared_ptr<LRadRt> lrad_ps(new LRadRt(rt_lrad_only->stokes_coefficient(),
                                                        config_atmosphere,
                                                        config_spectral_window->spectral_bound(),
                                                        sza, zen_small, azm, pure_nadir,
                                                        rt_lrad_only->number_stokes(),
                                                        false, 4, 0.01,
                                                        LRadDriver::REGULAR));
  std::string d = serialize_write_string(lrad_ps);
  if(false)
    std::cerr << d;
  boost::shared_ptr<LRadRt> lrad_ps_r =
    serialize_read_string<LRadRt>(d);
  Array<double, 1> wn(1);
  wn = 12929.94;
  Array<double, 1> stokes_expect(3);
  stokes_expect = 
    0.0088645271123479803949, -0.0011307484227398416094, -0.00041155876643655497481;
  BOOST_CHECK_MATRIX_CLOSE_TOL(lrad_ps_r->stokes(wn, 0)(0, Range::all()), stokes_expect, 1e-6);
}

BOOST_AUTO_TEST_CASE(plane_parallel_ps)
{
  // Set up to run l_rad in regular ps mode
  Array<double, 1> zen_small(3);
  zen_small = 1.0e-6;
  boost::shared_ptr<LRadRt> lrad_ps(new LRadRt(rt_lrad_only->stokes_coefficient(),
                                                        config_atmosphere,
                                                        config_spectral_window->spectral_bound(),
                                                        sza, zen_small, azm, pure_nadir,
                                                        rt_lrad_only->number_stokes(),
                                                        false, 4, 0.01,
                                                        LRadDriver::PLANE_PARALLEL));
  Array<double, 1> wn(1);
  wn = 12929.94;
  Array<double, 1> stokes_expect(3);
  stokes_expect = 
    0.0088439097970879985283, -0.0011295516868515330204, -0.0004111231901969267309;
  BOOST_CHECK_MATRIX_CLOSE_TOL(lrad_ps->stokes(wn, 0)(0, Range::all()), stokes_expect, 1e-6);
}

BOOST_AUTO_TEST_CASE(plane_parallel_ps_serialization)
{
  if(!have_serialize_supported())
    return;
  // Set up to run l_rad in regular ps mode
  Array<double, 1> zen_small(3);
  zen_small = 1.0e-6;
  boost::shared_ptr<LRadRt> lrad_ps(new LRadRt(rt_lrad_only->stokes_coefficient(),
                                                        config_atmosphere,
                                                        config_spectral_window->spectral_bound(),
                                                        sza, zen_small, azm, pure_nadir,
                                                        rt_lrad_only->number_stokes(),
                                                        false, 4, 0.01,
                                                        LRadDriver::PLANE_PARALLEL));
  std::string d = serialize_write_string(lrad_ps);
  if(false)
    std::cerr << d;
  boost::shared_ptr<LRadRt> lrad_ps_r =
    serialize_read_string<LRadRt>(d);
  Array<double, 1> wn(1);
  wn = 12929.94;
  Array<double, 1> stokes_expect(3);
  stokes_expect = 
    0.0088439097970879985283, -0.0011295516868515330204, -0.0004111231901969267309;
  BOOST_CHECK_MATRIX_CLOSE_TOL(lrad_ps_r->stokes(wn, 0)(0, Range::all()), stokes_expect, 1e-6);
}

BOOST_AUTO_TEST_SUITE_END()

// --------- coxmunk ----------- //

BOOST_FIXTURE_TEST_SUITE(l_rad_rt_coxmunk, LRadCoxmunkFixture)

BOOST_AUTO_TEST_CASE(reflectance_first_order)
{
  Array<double, 1> wn(2);
  wn = 12929.94, 12930.30;
  Array<double, 2> stokes_expect(2, 3);
  stokes_expect = 
    0.0049435008885970674331, -0.0013395752071591116635, -0.00048224005832755750985,
    0.0049438145263468770418, -0.0013397023526495689585, -0.00048228583003817615277;
  BOOST_CHECK_MATRIX_CLOSE_TOL(rt_first_order->stokes(wn, 0), 
                               stokes_expect, check_tol);
  BOOST_CHECK_MATRIX_CLOSE_TOL
    (rt_first_order->stokes_and_jacobian(wn, 0).value(),
     stokes_expect, check_tol);
  Array<double, 1> rad_expect(stokes_expect(Range::all(), 0));
  BOOST_CHECK_MATRIX_CLOSE_TOL(rt_first_order->reflectance(wn, 0, true).
                               spectral_range().data(), 
                               rad_expect, check_tol);
  BOOST_CHECK_MATRIX_CLOSE_TOL
    (rt_first_order->reflectance(wn, 0).spectral_range().data(), 
     rad_expect, check_tol);
}

BOOST_AUTO_TEST_CASE(reflectance_first_order_serialization)
{
  if(!have_serialize_supported())
    return;
  std::string d = serialize_write_string(rt_first_order);
  if(false)
    std::cerr << d;
  boost::shared_ptr<LRadRt> rt_first_order_r =
    serialize_read_string<LRadRt>(d);
  Array<double, 1> wn(2);
  wn = 12929.94, 12930.30;
  Array<double, 2> stokes_expect(2, 3);
  stokes_expect = 
    0.0049435008885970674331, -0.0013395752071591116635, -0.00048224005832755750985,
    0.0049438145263468770418, -0.0013397023526495689585, -0.00048228583003817615277;
  BOOST_CHECK_MATRIX_CLOSE_TOL(rt_first_order_r->stokes(wn, 0), 
                               stokes_expect, check_tol);
  BOOST_CHECK_MATRIX_CLOSE_TOL
    (rt_first_order_r->stokes_and_jacobian(wn, 0).value(),
     stokes_expect, check_tol);
  Array<double, 1> rad_expect(stokes_expect(Range::all(), 0));
  BOOST_CHECK_MATRIX_CLOSE_TOL(rt_first_order_r->reflectance(wn, 0, true).
                               spectral_range().data(), 
                               rad_expect, check_tol);
  BOOST_CHECK_MATRIX_CLOSE_TOL
    (rt_first_order_r->reflectance(wn, 0).spectral_range().data(), 
     rad_expect, check_tol);
}

BOOST_AUTO_TEST_CASE(reflectance_second_order)
{
  Array<double, 1> wn(2);
  wn = 12929.94, 12930.30;
  Array<double, 2> stokes_expect(2, 3);
  stokes_expect = 
    0.0049757072246467068388, -0.0017906853636642438773, -0.00061963297356329384482,
    0.0049760200639414903928, -0.0017908307512102411303, -0.0006196862328839785293;
   BOOST_CHECK_MATRIX_CLOSE_TOL(rt_second_order->stokes(wn, 0), 
                               stokes_expect, check_tol);
  BOOST_CHECK_MATRIX_CLOSE_TOL
    (rt_second_order->stokes_and_jacobian(wn, 0).value(), 
     stokes_expect, check_tol);
  Array<double, 1> rad_expect(stokes_expect(Range::all(), 0));
  BOOST_CHECK_MATRIX_CLOSE_TOL(rt_second_order->reflectance(wn, 0, true).
                               spectral_range().data(), 
                               rad_expect, check_tol);
  BOOST_CHECK_MATRIX_CLOSE_TOL
    (rt_second_order->reflectance(wn, 0).spectral_range().data(), 
     rad_expect, check_tol);
}

BOOST_AUTO_TEST_CASE(reflectance_second_order_serialization)
{
  if(!have_serialize_supported())
    return;
  std::string d = serialize_write_string(rt_second_order);
  if(false)
    std::cerr << d;
  boost::shared_ptr<LRadRt> rt_second_order_r =
    serialize_read_string<LRadRt>(d);
  Array<double, 1> wn(2);
  wn = 12929.94, 12930.30;
  Array<double, 2> stokes_expect(2, 3);
  stokes_expect = 
    0.0049757072246467068388, -0.0017906853636642438773, -0.00061963297356329384482,
    0.0049760200639414903928, -0.0017908307512102411303, -0.0006196862328839785293;
  BOOST_CHECK_MATRIX_CLOSE_TOL(rt_second_order_r->stokes(wn, 0), 
                               stokes_expect, check_tol);
  BOOST_CHECK_MATRIX_CLOSE_TOL
    (rt_second_order_r->stokes_and_jacobian(wn, 0).value(), 
     stokes_expect, check_tol);
  Array<double, 1> rad_expect(stokes_expect(Range::all(), 0));
  BOOST_CHECK_MATRIX_CLOSE_TOL(rt_second_order_r->reflectance(wn, 0, true).
                               spectral_range().data(), 
                               rad_expect, check_tol);
  BOOST_CHECK_MATRIX_CLOSE_TOL
    (rt_second_order_r->reflectance(wn, 0).spectral_range().data(), 
     rad_expect, check_tol);
}

BOOST_AUTO_TEST_CASE(jac_first_order)
{
  // We use rt_lrad_only so we only test the portion of the Jacobian coming
  // from LRad. We separately test LIDORT.
  is_long_test();
  StateVector& sv = *config_state_vector;
  Array<double, 1> sv0(config_initial_guess->initial_guess());
  sv.update_state(sv0);

  int spec_index = 2;
  Array<double, 1> wn(1);
  wn = 4820.0;

  // Test jacobians for all three l_rad_first pseudo spherical modes
  for(int ps_idx = LRadDriver::REGULAR; ps_idx <= LRadDriver::PLANE_PARALLEL; ps_idx++) {
    LRadDriver::PsMode ps_mode = static_cast<LRadDriver::PsMode>(ps_idx);
    boost::shared_ptr<LRadRt> lrad_ps(new LRadRt(rt_lrad_only->stokes_coefficient(),
                                                         config_atmosphere,
                                                         config_spectral_window->spectral_bound(),
                                                         sza, zen, azm, pure_nadir,
                                                         rt_lrad_only->number_stokes(),
                                                         false, 4, 0.01,
                                                         ps_mode));

    ArrayAd<double, 2> stk = 
      lrad_ps->stokes_and_jacobian(wn, spec_index);
    Array<double, 2> stk0(stk.shape());
    stk0 = stk.value();
    Array<double, 3> jac = stk.jacobian().copy();
    for(int i = 103; i < sv.state().rows(); ++i) {
    //  for(int i = 0; i < sv.state().rows(); ++i) {
      Array<double, 1> svn(sv0.copy());
      svn(i) += epsilon(i);
      sv.update_state(svn);
      Array<double, 2> jacfd(stk.shape());
      jacfd = (lrad_ps->stokes(wn, spec_index) - stk0) 
        / epsilon(i);
      double diff = max(abs(jac(Range::all(), Range::all(), i) - jacfd));
      // Can turn this on to dump values,
      // if needed for debugging
      bool debug = false;
      if(debug || (diff > 1e-6)) {
        std::cerr << "SV:" << i << ", Initial guess: " << sv0(i) << std::endl
                  << "epsilon = " << epsilon(i) << std::endl
                  << "abs diff = " << diff << std::endl;
        double fdmax_diff = max(abs(jacfd));
        if (fdmax_diff > 0.0)
          std::cerr << "rel diff = "<< diff / fdmax_diff << std::endl;
        std::cerr << "analytic: " << jac(0,Range::all(), i) << std::endl
                  << "finite diff: " << jacfd << std::endl;
      }
      BOOST_CHECK(diff < 1e-6);
    }
  }
}

BOOST_AUTO_TEST_CASE(jac_second_order)
{
  // We use rt_lrad_only so we only test the portion of the Jacobian coming
  // from LRad. We separately test LIDORT.
  is_long_test();
  StateVector& sv = *config_state_vector;
  Array<double, 1> sv0(config_initial_guess->initial_guess());
  sv.update_state(sv0);
  int spec_index = 2;
  Array<double, 1> wn(1);
  wn = 4820.0;
  ArrayAd<double, 2> stk = 
    rt_lrad_only_second->stokes_and_jacobian(wn, spec_index);
  Array<double, 2> stk0(stk.shape());
  stk0 = stk.value();
  Array<double, 3> jac = stk.jacobian().copy();
  for(int i = 0; i < sv.state().rows(); ++i) {
    Array<double, 1> svn(sv0.copy());
    svn(i) += epsilon(i);
    sv.update_state(svn);
    Array<double, 2> jacfd(stk.shape());
    jacfd = (rt_lrad_only_second->stokes(wn, spec_index) - stk0) 
      / epsilon(i);
    double diff = max(abs(jac(Range::all(), Range::all(), i) - jacfd));
    // Can turn this on to dump values,
    // if needed for debugging
    bool debug = false;
    if(debug || (diff > 1e-4)) {
      std::cerr << "SV:" << i << ", Initial guess: " << sv0(i) << std::endl
                << "epsilon = " << epsilon(i) << std::endl
                << "abs diff = " << diff << std::endl;
      double fdmax_diff = max(abs(jacfd));
      if (fdmax_diff > 0.0)
        std::cerr << "rel diff = "<< diff / fdmax_diff << std::endl;
      std::cerr << "analytic: " << jac(0,Range::all(), i) << std::endl
                << "finite diff: " << jacfd << std::endl;
    }
    BOOST_CHECK(diff < 1e-4);
  }
}

BOOST_AUTO_TEST_CASE(regular_ps)
{
  // Set up to run l_rad in regular ps mode
  Array<double, 1> zen_small(3);
  zen_small = 1.0e-6;
  boost::shared_ptr<LRadRt> lrad_ps(new LRadRt(rt_lrad_only->stokes_coefficient(),
                                                        config_atmosphere,
                                                        config_spectral_window->spectral_bound(),
                                                        sza, zen_small, azm, pure_nadir,
                                                        rt_lrad_only->number_stokes(),
                                                        false, 4, 0.01,
                                                        LRadDriver::REGULAR));
  Array<double, 1> wn(1);
  wn = 12929.94;
  Array<double, 1> stokes_expect(3);
  stokes_expect = 
    0.001853123287241498305, -0.0011307484227534182059, -0.00041155876643241315989;
  BOOST_CHECK_MATRIX_CLOSE_TOL(lrad_ps->stokes(wn, 0)(0, Range::all()), stokes_expect, 1e-6);
}

BOOST_AUTO_TEST_CASE(regular_ps_serialization)
{
  if(!have_serialize_supported())
    return;
  // Set up to run l_rad in regular ps mode
  Array<double, 1> zen_small(3);
  zen_small = 1.0e-6;
  boost::shared_ptr<LRadRt> lrad_ps(new LRadRt(rt_lrad_only->stokes_coefficient(),
                                                        config_atmosphere,
                                                        config_spectral_window->spectral_bound(),
                                                        sza, zen_small, azm, pure_nadir,
                                                        rt_lrad_only->number_stokes(),
                                                        false, 4, 0.01,
                                                        LRadDriver::REGULAR));
  std::string d = serialize_write_string(lrad_ps);
  if(false)
    std::cerr << d;
  boost::shared_ptr<LRadRt> lrad_ps_r =
    serialize_read_string<LRadRt>(d);
  Array<double, 1> wn(1);
  wn = 12929.94;
  Array<double, 1> stokes_expect(3);
  stokes_expect = 
    0.001853123287241498305, -0.0011307484227534182059, -0.00041155876643241315989;
  BOOST_CHECK_MATRIX_CLOSE_TOL(lrad_ps_r->stokes(wn, 0)(0, Range::all()), stokes_expect, 1e-6);
}

BOOST_AUTO_TEST_CASE(plane_parallel_ps)
{
  // Set up to run l_rad in regular ps mode
  Array<double, 1> zen_small(3);
  zen_small = 1.0e-6;
  boost::shared_ptr<LRadRt> lrad_ps(new LRadRt(rt_lrad_only->stokes_coefficient(),
                                                        config_atmosphere,
                                                        config_spectral_window->spectral_bound(),
                                                        sza, zen_small, azm, pure_nadir, 
                                                        rt_lrad_only->number_stokes(),
                                                        false, 4, 0.01,
                                                        LRadDriver::PLANE_PARALLEL));
  Array<double, 1> wn(1);
  wn = 12929.94;
  Array<double, 1> stokes_expect(3);
  stokes_expect = 
    0.0018510852209117484562, -0.0011295516868650760066, -0.00041112319019279640851;
  BOOST_CHECK_MATRIX_CLOSE_TOL(lrad_ps->stokes(wn, 0)(0, Range::all()), stokes_expect, 1e-6);
}

BOOST_AUTO_TEST_CASE(plane_parallel_ps_serialization)
{
  if(!have_serialize_supported())
    return;
  // Set up to run l_rad in regular ps mode
  Array<double, 1> zen_small(3);
  zen_small = 1.0e-6;
  boost::shared_ptr<LRadRt> lrad_ps(new LRadRt(rt_lrad_only->stokes_coefficient(),
                                                        config_atmosphere,
                                                        config_spectral_window->spectral_bound(),
                                                        sza, zen_small, azm, pure_nadir, 
                                                        rt_lrad_only->number_stokes(),
                                                        false, 4, 0.01,
                                                        LRadDriver::PLANE_PARALLEL));
  std::string d = serialize_write_string(lrad_ps);
  if(false)
    std::cerr << d;
  boost::shared_ptr<LRadRt> lrad_ps_r =
    serialize_read_string<LRadRt>(d);
  Array<double, 1> wn(1);
  wn = 12929.94;
  Array<double, 1> stokes_expect(3);
  stokes_expect = 
    0.0018510852209117484562, -0.0011295516868650760066, -0.00041112319019279640851;
  BOOST_CHECK_MATRIX_CLOSE_TOL(lrad_ps_r->stokes(wn, 0)(0, Range::all()), stokes_expect, 1e-6);
}

BOOST_AUTO_TEST_SUITE_END()

#include "lsi_rt.h"
#include "fp_logger.h"
#include "lidort_fixture.h"
#include "unit_test_support.h"
#include "hdf_file.h"

using namespace FullPhysics;
using namespace blitz;
BOOST_FIXTURE_TEST_SUITE(lsi_rt, LidortLowHighLambertianFixture)

BOOST_AUTO_TEST_CASE(err_est)
{
  is_long_test();                // Skip unless we are running long tests.
  turn_on_logger();                // Have log output show up.
  
  HdfFile config(test_data_dir() + "lua/example_static_input.h5");
  IfstreamCs expected_data(test_data_dir() + "expected/lsi_rt/lsi_expected");
  LsiRt rt(low_rt, high_rt, config);
  int wn_i = 0;
  for(double wn = 12929.94; wn <= 13210.15; wn += 0.01)
    wn_i += 1;
  Array<double, 1> wn_arr(wn_i);
  wn_i = 0;
  for(double wn = 12929.94; wn <= 13210.15; wn += 0.01, ++wn_i)
    wn_arr(wn_i) = wn;
  ArrayAd<double, 2> err_est = rt.correction_only(wn_arr, 0);
  Logger::info() << low_rt->atmosphere_ptr()->timer_info();
  Array<double, 1> err_est_expect;
  expected_data >> err_est_expect;
  if(false) {                   // Write to output, if we need to
                                // regenerate expected values.
    ofstream out_expt(test_data_dir() + "expected/lsi_rt/lsi_expected");
    out_expt << std::setprecision(20) << std::scientific 
             << err_est.value()(Range::all(), 0) << std::endl;
  }
  if(false) {
    for(int i = 0; i < err_est_expect.rows(); i++) {
      std::cerr << "[" << i << "]: e: " << err_est_expect(i) << " c: " << err_est.value()(i, 0) 
                << " d: " << (err_est_expect(i) - err_est.value()(i, 0)) << std::endl;
    }
  }
  BOOST_CHECK_MATRIX_CLOSE_TOL(err_est.value()(Range::all(), 0), 
                               err_est_expect, 1e-7);
}

BOOST_AUTO_TEST_CASE(err_est_jac)
{
  // Even for a long test, this takes a long time to run (several
  // minutes). We don't normally run this, since all the jacobian
  // calculation are done by AutoDerivative anyways which
  // automatically calculates this. But we can turn this test on if
  // needed.
  is_really_long_test();     // Skip unless we are running long tests.

  HdfFile c(test_data_dir() + "lua/example_static_input.h5");
  IfstreamCs expected_data(test_data_dir() + "expected/lsi_rt/lsi_expected");
  LsiRt rt(low_rt, high_rt, c);

  StateVector& sv = *config_state_vector;
  Array<double, 1> sv0(config_initial_guess->initial_guess());
  sv.update_state(sv0);

  int spec_index = 0;
  int wn_i = 0;
  for(double wn = 13000.0; wn <= 13000.1; wn += 0.01)
    wn_i += 1;
  Array<double, 1> wn_arr(wn_i);
  wn_i = 0;
  for(double wn = 13000.0; wn <= 13000.1; wn += 0.01, ++wn_i)
    wn_arr(wn_i) = wn;

  ArrayAd<double, 2> err_est = rt.correction_only(wn_arr, spec_index);
  Array<double, 2> e0(err_est.shape());
  e0 = err_est.value();
  Array<double, 3> jac = err_est.jacobian().copy();
  for(int i = 0; i < sv.state().rows(); ++i) {
    Array<double, 1> svn(sv0.copy());
    svn(i) += epsilon(i);
    sv.update_state(svn);
    Array<double, 2> jacfd(err_est.shape());
    jacfd = (rt.correction_only(wn_arr, spec_index).value() - e0) / epsilon(i);
    Array<double, 2> diff(jac(Range::all(), Range::all(), i) - jacfd);
    if(false) {                 // Can turn this on to dump values,
                                // if needed for debugging
      if(max(abs(diff)) !=0) {
        std::cerr << i << ": " 
            << "jac_a = " << jac(Range::all(), Range::all(), i) << std::endl
            << "jac_fd = " << jacfd << std::endl
            << "abs_diff = " << max(abs(diff)) << " "
            << "rel_diff = " << max(abs(where(abs(jacfd) < 1e-15, 0, diff/jacfd))) << "\n";
      }
    }
    BOOST_CHECK(max(abs(diff)) < 3e-4);
  }
}

BOOST_AUTO_TEST_CASE(stokes)
{
  // We don't normally need to run this test. All the real calculation
  // is in err_est test above. Stokes is just a direct function of the
  // err_est. We need a test to make sure that the actual calculation
  // is ok, but given the time it takes to run this we don't need to
  // generally do so.
  is_really_long_test();     // Skip unless we are running long tests.
  turn_on_logger();                // Have log output show up.
  
  HdfFile config(test_data_dir() + "lua/example_static_input.h5");
  IfstreamCs expected_data(test_data_dir() + "expected/lsi_rt/stokes");
  LsiRt rt(low_rt, high_rt, config);
  int wn_i = 0;
  for(double wn = 12929.94; wn <= 13210.15; wn += 0.01)
    wn_i += 1;
  Array<double, 1> wn_arr(wn_i);
  wn_i = 0;
  for(double wn = 12929.94; wn <= 13210.15; wn += 0.01, ++wn_i)
    wn_arr(wn_i) = wn;
  Array<double, 1> refl = rt.reflectance(wn_arr, 0, true).spectral_range().data();
  Array<double, 1> refl_expect;
  expected_data >> refl_expect;
  BOOST_CHECK_MATRIX_CLOSE(refl(Range(9,19)), refl_expect);
  expected_data >> refl_expect;
  BOOST_CHECK_MATRIX_CLOSE(refl(Range(15999,16009)), refl_expect);
  if(false) {                   // Write to output, if we need to
                                // regenerate expected values.
    std::ofstream out_expt(test_data_dir() + "expected/lsi_rt/stokes");
    out_expt << std::setprecision(20) << std::scientific
             << "# refl Range(9,19)" << std::endl
             << std::scientific << refl(Range(9,19)) << std::endl
             << "# refl Range(15999,16009)" << std::endl
             << std::scientific << refl(Range(15999,16009)) << std::endl;
  }

}

BOOST_AUTO_TEST_CASE(stokes_and_jacobian)
{
  // We don't normally need to run this test. All the real calculation
  // is in err_est test above. Stokes is just a direct function of the
  // err_est. We need a test to make sure that the actual calculation
  // is ok, but given the time it takes to run this we don't need to
  // generally do so.
  is_really_long_test();     // Skip unless we are running long tests.
  turn_on_logger();                // Have log output show up.

  HdfFile config(test_data_dir() + "lua/example_static_input.h5");
  IfstreamCs expected_data(test_data_dir() + "expected/lsi_rt/stokes");
  LsiRt rt(low_rt, high_rt, config);
  int wn_i = 0;
  for(double wn = 12929.94; wn <= 13210.15; wn += 0.01)
    wn_i += 1;
  Array<double, 1> wn_arr(wn_i);
  wn_i = 0;
  for(double wn = 12929.94; wn <= 13210.15; wn += 0.01, ++wn_i)
    wn_arr(wn_i) = wn;
  ArrayAd<double, 1> refl = rt.reflectance(wn_arr, 0).spectral_range().data_ad();
  Array<double, 1> refl_expect;
  expected_data >> refl_expect;
  BOOST_CHECK_MATRIX_CLOSE(refl.value()(Range(9,19)), refl_expect);
  expected_data >> refl_expect;
  BOOST_CHECK_MATRIX_CLOSE(refl.value()(Range(15999,16009)), refl_expect);
  if(false) {                        // Write to output, if we need to
                                // regenerate expected values.
    std::cerr.precision(20);
    std::cerr << "# refl Range(9,19)" << std::endl
              << std::scientific << refl.value()(Range(9,19)) << std::endl
              << "# refl Range(15999,16009)" << std::endl
              << std::scientific << refl.value()(Range(15999,16009)) << std::endl;
  }

}

// This next text is in response to Ticket #939. When we run the LSI
// twice on the same input, we expect to get the same results. For the
// very first point (and only for the very first point), we had a
// condition where the Jacobians were difference, primarily for the
// CO2 VMR. This test illustrates this problem.
BOOST_AUTO_TEST_CASE(lsi_run_twice)
{
  is_long_test();              // Skip unless we are running long tests.
  turn_on_logger();            // Have log output show up.
  int band = 1;                       // Look at weak CO2 
  SpectralDomain spec_domain = highres_grid(band);
  ArrayAd<double, 1> refl = config_rt->reflectance(spec_domain,band).spectral_range().data_ad();
  ArrayAd<double, 1> refl2 = config_rt->reflectance(spec_domain,band).spectral_range().data_ad();
  BOOST_CHECK_MATRIX_CLOSE(refl.jacobian(), refl2.jacobian());
}

BOOST_AUTO_TEST_SUITE_END()

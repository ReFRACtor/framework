#include "error_analysis.h"
#include "absorber_absco.h"
#include "unit_test_support.h"
#include "solver_finished_fixture.h"

using namespace FullPhysics;
using namespace blitz;
BOOST_FIXTURE_TEST_SUITE(error_analysis, SolverFinishedFixture)

BOOST_AUTO_TEST_CASE(basic)
{
  is_long_test();                // Skip unless we are running long tests.
  IfstreamCs expected_data(test_data_dir() + "expected/error_analysis/basic");
  const ErrorAnalysis& err = *config_error_analysis;
  BOOST_CHECK_CLOSE(err.signal_level(0), 6.3001985485239222e+19, 1e-2);
  BOOST_CHECK_CLOSE(err.signal_level(1), 2.4251725004008743e+19, 1e-2);
  BOOST_CHECK_CLOSE(err.signal_level(2), 1.0586716647420277e+19, 1e-2);
  BOOST_CHECK_CLOSE(err.noise_level(0), 1.9566716780094205e+17, 1e-2);
  BOOST_CHECK_CLOSE(err.noise_level(1), 52415731785312176, 1e-2);
  BOOST_CHECK_CLOSE(err.noise_level(2), 37665172560823952, 1e-2);
  BOOST_CHECK_CLOSE(err.residual_sum_sq(0), 1.1369857341876807e-16, 1e-2);
  BOOST_CHECK_CLOSE(err.residual_sum_sq(1), 1.3128118455044732e-16, 1e-2);
  BOOST_CHECK_CLOSE(err.residual_sum_sq(2), 1.0780773740167574e-17, 1e-2);
  BOOST_CHECK_CLOSE(err.residual_mean_sq(0), 3.690067156436877e-10, 1e-2);
  BOOST_CHECK_CLOSE(err.residual_mean_sq(1), 4.4872208553479776e-10, 1e-2);
  BOOST_CHECK_CLOSE(err.residual_mean_sq(2), 1.152960791592836e-10, 1e-2);
  BOOST_CHECK_CLOSE(err.relative_residual_mean_sq(0), err.residual_mean_sq(0) / err.signal_level(0), 1e-2);
  BOOST_CHECK_CLOSE(err.relative_residual_mean_sq(1), err.residual_mean_sq(1) / err.signal_level(1), 1e-2);
  BOOST_CHECK_CLOSE(err.relative_residual_mean_sq(2), err.residual_mean_sq(2) / err.signal_level(2), 1e-2);
  BOOST_CHECK_CLOSE(err.reduced_chisq(0), 0.037459152030236656, 1e-2);
  BOOST_CHECK_CLOSE(err.reduced_chisq(2), 0.031025673876127577, 1e-2);
  BOOST_CHECK_CLOSE(err.reduced_chisq(1), 0.053007212109111294, 1e-2);
  BOOST_CHECK_CLOSE(err.degrees_of_freedom(), 15.762273734062212, 1e-4);
}
BOOST_AUTO_TEST_SUITE_END()

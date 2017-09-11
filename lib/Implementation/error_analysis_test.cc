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
  BOOST_CHECK_CLOSE(err.signal(0), 6.3001985485239222e+19, 1e-2);
  BOOST_CHECK_CLOSE(err.signal(1), 2.4251725004008743e+19, 1e-2);
  BOOST_CHECK_CLOSE(err.signal(2), 1.0586716647420277e+19, 1e-2);
  BOOST_CHECK_CLOSE(err.noise(0), 1.9566716780094205e+17, 1e-2);
  BOOST_CHECK_CLOSE(err.noise(1), 52415731785312176, 1e-2);
  BOOST_CHECK_CLOSE(err.noise(2), 37665172560823952, 1e-2);
  BOOST_CHECK_CLOSE(err.residual_sum_sq(0), 1.1369857341876807e-16, 1e-2);
  BOOST_CHECK_CLOSE(err.residual_sum_sq(1), 1.3128118455044732e-16, 1e-2);
  BOOST_CHECK_CLOSE(err.residual_sum_sq(2), 1.0780773740167574e-17, 1e-2);
  BOOST_CHECK_CLOSE(err.residual_mean_sq(0), 3.690067156436877e-10, 1e-2);
  BOOST_CHECK_CLOSE(err.residual_mean_sq(1), 4.4872208553479776e-10, 1e-2);
  BOOST_CHECK_CLOSE(err.residual_mean_sq(2), 1.152960791592836e-10, 1e-2);
  BOOST_CHECK_CLOSE(err.relative_residual_mean_sq(0), err.residual_mean_sq(0) / err.signal(0), 1e-2);
  BOOST_CHECK_CLOSE(err.relative_residual_mean_sq(1), err.residual_mean_sq(1) / err.signal(1), 1e-2);
  BOOST_CHECK_CLOSE(err.relative_residual_mean_sq(2), err.residual_mean_sq(2) / err.signal(2), 1e-2);
  BOOST_CHECK_CLOSE(err.reduced_chisq(0), 0.037459152030236656, 1e-2);
  BOOST_CHECK_CLOSE(err.reduced_chisq(2), 0.031025673876127577, 1e-2);
  BOOST_CHECK_CLOSE(err.reduced_chisq(1), 0.053007212109111294, 1e-2);
  BOOST_CHECK_CLOSE(err.xco2_uncert_noise(), sqrt(err.xco2_measurement_error()), 1e-4);
  BOOST_CHECK_CLOSE(err.xco2_uncert_smooth(), sqrt(err.xco2_smoothing_error()), 1e-4);
  BOOST_CHECK_CLOSE(err.xco2_uncert_interf(), sqrt(err.xco2_interference_error()), 1e-4);
  BOOST_CHECK_CLOSE(err.xco2_measurement_error(), 1.130188545667912e-11, 1e-4);
  BOOST_CHECK_CLOSE(err.xco2_smoothing_error(), 1.0785813873655232e-11, 1e-4);
  BOOST_CHECK_CLOSE(err.xco2_interference_error(), 5.9547559701327288e-12, 1e-4);
  // Equation 3-109 in the ATB
  BOOST_CHECK_CLOSE(err.xco2_uncertainty(), 4.6652613837991292e-06, 1e-4);
  BOOST_CHECK_CLOSE(err.degrees_of_freedom_full_vector(), 15.762273734062212, 1e-4);
  BOOST_CHECK_CLOSE(err.degrees_of_freedom_xco2(), 0.91615188124993685, 1e-4);

  Array<double, 1> xco2_avg_kernel_expected, xco2_avg_kernel_norm_expected,
    xco2_correlation_interf_expected, isu_expected, xco2_gain_expected;
  expected_data >> xco2_avg_kernel_expected
                >> xco2_avg_kernel_norm_expected
                >> xco2_correlation_interf_expected
                >> isu_expected
                >> xco2_gain_expected;
  BOOST_CHECK_MATRIX_CLOSE(err.xco2_avg_kernel(), xco2_avg_kernel_expected);
  BOOST_CHECK_MATRIX_CLOSE_TOL(err.xco2_avg_kernel_norm(), 
                               xco2_avg_kernel_norm_expected, 1e-6);
  BOOST_CHECK_MATRIX_CLOSE(err.xco2_correlation_interf(), 
                           xco2_correlation_interf_expected);
  BOOST_CHECK_MATRIX_CLOSE_TOL(err.interference_smoothing_uncertainty(),
                               isu_expected, 1e-13);
  BOOST_CHECK_MATRIX_CLOSE_TOL(err.xco2_gain_vector(), xco2_gain_expected, 1e-8);

}
BOOST_AUTO_TEST_SUITE_END()

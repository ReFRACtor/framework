#include "unit_test_support.h"
#include "error_analysis_output.h"
#include "output_hdf.h"
#include "solver_finished_fixture.h"

using namespace FullPhysics;
using namespace blitz;
BOOST_FIXTURE_TEST_SUITE(error_analysis_output, SolverFinishedFixture)

BOOST_AUTO_TEST_CASE(basic)
{
  is_long_test();		// Skip unless we are running long tests.

  blitz::Array<bool, 1> spec_flag(config_forward_model->num_channels());
  spec_flag = true;

  std::vector<std::string> band_name;
  band_name.push_back("A-Band");
  band_name.push_back("Weak-CO2");
  band_name.push_back("Strong-CO2");

  ErrorAnalysisOutput po(config_error_analysis, spec_flag, band_name, true);
  boost::shared_ptr<OutputHdf> out(new OutputHdf("error_analysis_output.h5", 20, 112, 5, 3));
  add_file_to_cleanup("error_analysis_output.h5");
  po.register_output(out);

  // Simple test, we just make sure that we can write output. All the
  // actual value calculation is checked in ErrorAnalysis unit test.

  out->write();
  BOOST_CHECK(true);
}

BOOST_AUTO_TEST_SUITE_END()



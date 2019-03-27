#include "configuration_fixture.h"
#include "unit_test_support.h"
#include "observation_level_1b.h"

using namespace FullPhysics;
using namespace blitz;

BOOST_FIXTURE_TEST_SUITE(observation_level_1b, ConfigurationFixture)

BOOST_AUTO_TEST_CASE(basic)
{
    boost::shared_ptr<ForwardModelSpectralGrid> fg(new ForwardModelSpectralGrid(config_instrument, config_spectral_window, config_spectrum_sampling));

    ObservationLevel1b obs_l1b(config_level_1b, config_instrument, fg);

    Spectrum radiance_all = obs_l1b.radiance_all();

    IfstreamCs expected_input(test_data_dir() + "expected/observation_level_1b/radiance_all");

    Array<double, 1> expected_domain, expected_range;
    expected_input >> expected_domain >> expected_range;
    if(false) {
      std::cerr.precision(20);
      std::cerr << "# all - spectral_domain\n"
		<< radiance_all.spectral_domain().data()
		<< "# all - spectral_range\n"
		<< radiance_all.spectral_range().data();
    }
    BOOST_CHECK_MATRIX_CLOSE_TOL(radiance_all.spectral_domain().data(), expected_domain, 1e-10);
    BOOST_CHECK_MATRIX_CLOSE_TOL(radiance_all.spectral_range().data(), expected_range, 1e-1);
}

BOOST_AUTO_TEST_SUITE_END()

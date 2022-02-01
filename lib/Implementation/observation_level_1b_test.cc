#include "serialized_configuration_fixture.h"
#include "unit_test_support.h"
#include "observation_level_1b.h"

using namespace FullPhysics;
using namespace blitz;

BOOST_FIXTURE_TEST_SUITE(observation_level_1b, LambertianConfigurationFixture)

void check_obs_l1b(const boost::shared_ptr<ObservationLevel1b>& obs_l1b, const string& test_data_dir) {
    Spectrum radiance_all = obs_l1b->radiance_all();

    IfstreamCs expected_input(test_data_dir + "expected/observation_level_1b/radiance_all");

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

BOOST_AUTO_TEST_CASE(basic)
{
    boost::shared_ptr<ForwardModelSpectralGrid> fg(new ForwardModelSpectralGrid(config_instrument, config_spectral_window, config_spectrum_sampling));

    boost::shared_ptr<ObservationLevel1b> obs_l1b(new ObservationLevel1b(config_level_1b, config_instrument, fg));

    check_obs_l1b(obs_l1b, test_data_dir());
}

BOOST_AUTO_TEST_CASE(serialization)
{
  if(!have_serialize_supported())
    return;

  boost::shared_ptr<ForwardModelSpectralGrid> fg(new ForwardModelSpectralGrid(config_instrument, config_spectral_window, config_spectrum_sampling));

  boost::shared_ptr<ObservationLevel1b> obs_l1b_orig(new ObservationLevel1b(config_level_1b, config_instrument, fg));

  std::string serial_str = serialize_write_string(obs_l1b_orig);
  boost::shared_ptr<ObservationLevel1b> obs_l1b_read = serialize_read_string<ObservationLevel1b>(serial_str);

  check_obs_l1b(obs_l1b_read, test_data_dir());
}

BOOST_AUTO_TEST_SUITE_END()

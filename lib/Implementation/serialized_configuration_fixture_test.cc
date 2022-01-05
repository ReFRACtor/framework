#include "unit_test_support.h"

#include "serialized_configuration_fixture.h"

using namespace blitz;
using namespace FullPhysics;

BOOST_FIXTURE_TEST_SUITE(serialized_configuration_fixture, GlobalFixture)

BOOST_AUTO_TEST_CASE(read)
{
    if(!have_serialize_supported())
        return;

    SerializedConfigurationFixture fixture(test_data_dir() + "in/configuration_fixture/example_base_config.bin.gz");

    // Check that objects are not empty
    BOOST_CHECK(fixture.config_absorber);
    BOOST_CHECK(fixture.config_pressure);
    BOOST_CHECK(fixture.config_aerosol);
    BOOST_CHECK(fixture.config_atmosphere);
    BOOST_CHECK(fixture.config_state_vector);
    BOOST_CHECK(fixture.config_solver);
    BOOST_CHECK(fixture.config_spectral_window);
    BOOST_CHECK(fixture.config_initial_guess);
    BOOST_CHECK(fixture.config_instrument);
    BOOST_CHECK(fixture.config_temperature);
    BOOST_CHECK(fixture.config_spectrum_sampling);
    BOOST_CHECK(fixture.config_level_1b);
    BOOST_CHECK(fixture.config_ground);
    BOOST_CHECK(fixture.config_forward_model);
    BOOST_CHECK(fixture.config_observation);
    BOOST_CHECK(fixture.config_rt);

}

BOOST_AUTO_TEST_SUITE_END()

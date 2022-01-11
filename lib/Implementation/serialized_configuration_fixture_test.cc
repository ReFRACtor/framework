#include "unit_test_support.h"

#include "lua_configuration_fixture.h"
#include "serialized_configuration_fixture.h"

using namespace blitz;
using namespace FullPhysics;

int DEBUG = false;

BOOST_FIXTURE_TEST_SUITE(serialized_configuration_fixture, LambertianConfigurationFixture)

BOOST_AUTO_TEST_CASE(valid_objects)
{
    if(!have_serialize_supported())
        return;

    // Check that objects are not empty
    BOOST_CHECK(config_absorber);
    BOOST_CHECK(config_pressure);
    BOOST_CHECK(config_aerosol);
    BOOST_CHECK(config_atmosphere);
    BOOST_CHECK(config_state_vector);
    BOOST_CHECK(config_solver);
    BOOST_CHECK(config_spectral_window);
    BOOST_CHECK(config_initial_guess);
    BOOST_CHECK(config_instrument);
    BOOST_CHECK(config_temperature);
    BOOST_CHECK(config_spectrum_sampling);
    BOOST_CHECK(config_level_1b);
    BOOST_CHECK(config_ground);
    BOOST_CHECK(config_forward_model);
    BOOST_CHECK(config_observation);
    BOOST_CHECK(config_rt);

}

BOOST_AUTO_TEST_CASE(epsilon_check)
{
    blitz::Array<double, 1> epsilon_expected(config_state_vector->observer_claimed_size());

    // Legacy way of configuring epsilon with hard coded indexes
    epsilon_expected = 1e-6;                  // Default
    epsilon_expected(Range(0, 19)) = 1e-7;    // CO2 VMR
    epsilon_expected(21) = 1e-3;              // Surface Pressure
    epsilon_expected(22) = 1e-4;              // Temperature
    epsilon_expected(Range(23, 34)) = 1e-8;   // Ground + Dispersion

    // For debugging epsilon values
    if (DEBUG) {
        for (int sv_idx = 0 ; sv_idx < epsilon_expected.rows(); sv_idx++) {
            std::string sv_name = config_state_vector->state_vector_name()(sv_idx);

            std::cerr << sv_idx << ": " << sv_name << ", epsilon = " << epsilon(sv_idx) << ", epsilon_expected = " << epsilon_expected(sv_idx) << std::endl;
        }
    }

    BOOST_CHECK_MATRIX_CLOSE(epsilon_expected, epsilon);    
}

BOOST_AUTO_TEST_CASE(state_check)
{
    LuaConfigurationFixture lua_fixture = LuaConfigurationFixture();

    BOOST_CHECK_EQUAL(lua_fixture.config_state_vector->state().rows(), config_state_vector->state().rows());

    for(int sv_idx = 0; sv_idx < config_state_vector->state().rows(); sv_idx++) {
        if (DEBUG) {
            std::cerr << sv_idx << " Lua:    " << lua_fixture.config_state_vector->state_vector_name()(sv_idx) << " = "<< lua_fixture.config_state_vector->state()(sv_idx) << std::endl
                      << sv_idx << " Serial: " << config_state_vector->state_vector_name()(sv_idx) << " = " << config_state_vector->state()(sv_idx) << std::endl;
        }
        BOOST_CHECK_CLOSE(lua_fixture.config_state_vector->state()(sv_idx), config_state_vector->state()(sv_idx), 1e-8);
    }
}
 
BOOST_AUTO_TEST_SUITE_END()

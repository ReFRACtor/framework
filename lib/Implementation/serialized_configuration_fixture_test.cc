#include "unit_test_support.h"

#include "lua_configuration_fixture.h"
#include "serialized_configuration_fixture.h"

using namespace blitz;
using namespace FullPhysics;

int DEBUG = false;

void check_valid_objects(const ConfigurationFixture& fixture)
{
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

void check_epsilon(const ConfigurationFixture& fixture)
{
    blitz::Array<double, 1> epsilon_expected(fixture.config_state_vector->observer_claimed_size());

    // Legacy way of configuring epsilon with hard coded indexes
    epsilon_expected = 1e-6;                  // Default
    epsilon_expected(Range(0, 19)) = 1e-7;    // CO2 VMR
    epsilon_expected(21) = 1e-3;              // Surface Pressure
    epsilon_expected(22) = 1e-4;              // Temperature
    epsilon_expected(Range(23, 34)) = 1e-8;   // Ground + Dispersion

    // For debugging epsilon values
    if (DEBUG) {
        for (int sv_idx = 0 ; sv_idx < epsilon_expected.rows(); sv_idx++) {
            std::string sv_name = fixture.config_state_vector->state_vector_name()(sv_idx);

            std::cerr << sv_idx << ": " << sv_name << ", epsilon = " << fixture.epsilon(sv_idx) << ", epsilon_expected = " << epsilon_expected(sv_idx) << std::endl;
        }
    }

    BOOST_CHECK_MATRIX_CLOSE(epsilon_expected, fixture.epsilon);
}

void check_state(const SerializedConfigurationFixture& ser_fixture, const LuaConfigurationFixture& lua_fixture)
{
    BOOST_CHECK_EQUAL(lua_fixture.config_state_vector->state().rows(), ser_fixture.config_state_vector->state().rows());

    for(int sv_idx = 0; sv_idx < ser_fixture.config_state_vector->state().rows(); sv_idx++) {
        if (DEBUG) {
            std::cerr << sv_idx << " Lua:    " << lua_fixture.config_state_vector->state_vector_name()(sv_idx) << " = " << lua_fixture.config_state_vector->state()(sv_idx) << std::endl
                      << sv_idx << " Serial: " << ser_fixture.config_state_vector->state_vector_name()(sv_idx) << " = " << ser_fixture.config_state_vector->state()(sv_idx) << std::endl;
        }
        BOOST_CHECK_CLOSE(lua_fixture.config_state_vector->state()(sv_idx), ser_fixture.config_state_vector->state()(sv_idx), 1e-8);
    }
}

/****************************************************************//**
 Lambertian Fixture
*******************************************************************/

BOOST_FIXTURE_TEST_SUITE(lambertian_configuration_fixture, LambertianConfigurationFixture)

BOOST_AUTO_TEST_CASE(valid_objects)
{
    if(!have_serialize_supported())
        return;
    
    check_valid_objects(*this);
}

BOOST_AUTO_TEST_CASE(epsilon_check)
{
    if(!have_serialize_supported())
        return;

    check_epsilon(*this);
}

BOOST_AUTO_TEST_CASE(state_check)
{
    if(!have_serialize_supported())
        return;

    LuaConfigurationFixture lua_fixture = LuaConfigurationFixture();

    check_state(*this, lua_fixture);
}
 
BOOST_AUTO_TEST_SUITE_END()

/****************************************************************//**
 Coxmunk Fixture
*******************************************************************/

BOOST_FIXTURE_TEST_SUITE(coxmunk_configuration_fixture, CoxmunkConfigurationFixture)

BOOST_AUTO_TEST_CASE(valid_objects)
{
    if(!have_serialize_supported())
        return;
    
    check_valid_objects(*this);
}

BOOST_AUTO_TEST_CASE(epsilon_check)
{
    if(!have_serialize_supported())
        return;

    check_epsilon(*this);
}

BOOST_AUTO_TEST_CASE(state_check)
{
    if(!have_serialize_supported())
        return;

    LuaCoxmunkConfigurationFixture lua_fixture = LuaCoxmunkConfigurationFixture();

    check_state(*this, lua_fixture);
}
 
BOOST_AUTO_TEST_SUITE_END()

/****************************************************************//**
 Coxmunk Plus Lambertian Fixture
*******************************************************************/

BOOST_FIXTURE_TEST_SUITE(coxmunk_lambertian_configuration_fixture, CoxmunkPlusLambertianConfigurationFixture)

BOOST_AUTO_TEST_CASE(valid_objects)
{
    if(!have_serialize_supported())
        return;
    
    check_valid_objects(*this);
}

BOOST_AUTO_TEST_CASE(epsilon_check)
{
    if(!have_serialize_supported())
        return;

    check_epsilon(*this);
}

BOOST_AUTO_TEST_CASE(state_check)
{
    if(!have_serialize_supported())
        return;

    LuaCoxmunkPlusLambertianConfigurationFixture lua_fixture = LuaCoxmunkPlusLambertianConfigurationFixture();

    check_state(*this, lua_fixture);
}
 
BOOST_AUTO_TEST_SUITE_END()

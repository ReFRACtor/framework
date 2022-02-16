#include "unit_test_support.h"

#ifdef HAVE_LUA
#include "lua_configuration_fixture.h"
#endif

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

void compare_state_vector(const Array<double, 1>& sv_expt, const Array<double, 1>& sv_test, const Array<std::string, 1>& names_expt, const Array<std::string, 1>& names_test)
{
    BOOST_CHECK_EQUAL(sv_expt.rows(), sv_test.rows());

    for(int sv_idx = 0; sv_idx < sv_test.rows(); sv_idx++) {
        if (DEBUG) {
            std::cerr << sv_idx << " Lua:    " << names_expt(sv_idx) << " = " << sv_expt(sv_idx) << std::endl
                      << sv_idx << " Serial: " << names_test(sv_idx) << " = " << sv_test(sv_idx) << std::endl;
        }
        BOOST_CHECK_CLOSE(sv_expt(sv_idx), sv_test(sv_idx), 1e-8);
    }
}

#ifdef HAVE_LUA
void check_state_vs_lua(const SerializedConfigurationFixture& ser_fixture, const LuaConfigurationFixture& lua_fixture, const std::string& expected_name, const std::string& test_data_dir)
{
    
    // Overwrite expected values for when Lua is not enabled
    if (false) {
        ofstream out(test_data_dir + "expected/serialized_configuration_fixture/state_vector_" + expected_name);
        out << std::setprecision(20) << std::scientific;
        out << "# " << expected_name << " state_vector state" << std::endl;
        out << lua_fixture.config_state_vector->state() << std::endl;
        out << "# " << expected_name << " state_vector names" << std::endl;
        for(int sv_idx = 0; sv_idx < lua_fixture.config_state_vector->state().rows(); sv_idx++) {
            out << lua_fixture.config_state_vector->state_vector_name()(sv_idx) << std::endl;
        }
    }

    compare_state_vector(lua_fixture.config_state_vector->state(), ser_fixture.config_state_vector->state(), lua_fixture.config_state_vector->state_vector_name(), ser_fixture.config_state_vector->state_vector_name());
}
#endif

void check_state_vs_expected(const SerializedConfigurationFixture& ser_fixture, const std::string& expected_name, const std::string& test_data_dir)
{
    IfstreamCs expt_data(test_data_dir + "expected/serialized_configuration_fixture/state_vector_" + expected_name);
    
    Array<double, 1> state_expt;
    expt_data >> state_expt;

    // Required to transistion from stream operators to getline
    // Causes stream to ignore up until the next newline
    expt_data.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

    // Delete next empty line then comment line
    std::string discard;
    std::getline(expt_data, discard);
    std::getline(expt_data, discard);
       
    Array<std::string, 1> names_expt(state_expt.rows());
    for(int sv_idx = 0; sv_idx < state_expt.rows(); sv_idx++) {
        std::string line;
        std::getline(expt_data, line);
        names_expt(sv_idx) = line;
    }

    compare_state_vector(state_expt, ser_fixture.config_state_vector->state(), names_expt, ser_fixture.config_state_vector->state_vector_name());
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

#ifdef HAVE_LUA
    LuaConfigurationFixture lua_fixture = LuaConfigurationFixture();
    check_state_vs_lua(*this, lua_fixture, "lambertian", test_data_dir());
#endif
    check_state_vs_expected(*this, "lambertian", test_data_dir());
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

#ifdef HAVE_LUA
    LuaCoxmunkConfigurationFixture lua_fixture = LuaCoxmunkConfigurationFixture();
    check_state_vs_lua(*this, lua_fixture, "coxmunk", test_data_dir());
#endif
    check_state_vs_expected(*this, "coxmunk", test_data_dir());
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

#ifdef HAVE_LUA
    LuaCoxmunkPlusLambertianConfigurationFixture lua_fixture = LuaCoxmunkPlusLambertianConfigurationFixture();
    check_state_vs_lua(*this, lua_fixture, "coxmunk_lambertian", test_data_dir());
#endif
    check_state_vs_expected(*this, "coxmunk_lambertian", test_data_dir());
}

BOOST_AUTO_TEST_SUITE_END()

/****************************************************************//**
 BRDF Vegetation Fixture
*******************************************************************/

BOOST_FIXTURE_TEST_SUITE(brdf_veg_configuration_fixture, BrdfVegConfigurationFixture)

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

#ifdef HAVE_LUA
    LuaBrdfVegConfigurationFixture lua_fixture = LuaBrdfVegConfigurationFixture();
    check_state_vs_lua(*this, lua_fixture, "brdf_veg", test_data_dir());
#endif
    check_state_vs_expected(*this, "brdf_veg", test_data_dir());
}

BOOST_AUTO_TEST_SUITE_END()
 
/****************************************************************//**
 BRDF Soil Fixture
*******************************************************************/

BOOST_FIXTURE_TEST_SUITE(brdf_soil_configuration_fixture, BrdfSoilConfigurationFixture)

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

#ifdef HAVE_LUA
    LuaBrdfSoilConfigurationFixture lua_fixture = LuaBrdfSoilConfigurationFixture();
    check_state_vs_lua(*this, lua_fixture, "brdf_soil", test_data_dir());
#endif
    check_state_vs_expected(*this, "brdf_soil", test_data_dir());
}
 
BOOST_AUTO_TEST_SUITE_END()

/****************************************************************//**
 Two Broadener Fixture
*******************************************************************/

BOOST_FIXTURE_TEST_SUITE(two_broadener_configuration_fixture, TwoBroadenerConfigurationFixture)

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

#ifdef HAVE_LUA
    LuaTwoBroadenerConfigurationFixture lua_fixture = LuaTwoBroadenerConfigurationFixture();
    check_state_vs_lua(*this, lua_fixture, "two_broadener", test_data_dir());
#endif
    check_state_vs_expected(*this, "two_broadener", test_data_dir());
}
 
BOOST_AUTO_TEST_SUITE_END()

#include "lua_configuration_fixture.h"
#include <fenv.h>

using namespace FullPhysics;
using namespace blitz;

// Forward reference
std::map<std::string, boost::shared_ptr<LuaState> > LuaConfigurationFixture::config;

LuaConfigurationFixture::LuaConfigurationFixture(const std::string& Config_file)
  : config_filename(Config_file)
{
  init_variables();
  init_epsilon();
}

void LuaConfigurationFixture::init_variables()
{
  if(!config[config_filename]) {
    // Disable floating point exeptions while loading 
    // Lua configuration due to the way Lua parses certain
    // things.

    // Mac doesn't have this function, even though it is a C99
    // function. We check for this during configuration.
#ifdef HAVE_FEENABLEEXCEPT
    fedisableexcept(FE_INVALID | FE_DIVBYZERO | FE_OVERFLOW);
#endif
    config[config_filename] = LuaState::load_file(test_data_dir() + "/lua/" + config_filename);
#ifdef HAVE_FEENABLEEXCEPT
    feclearexcept(FE_INVALID | FE_DIVBYZERO | FE_OVERFLOW);
    feenableexcept(FE_INVALID | FE_DIVBYZERO | FE_OVERFLOW);
#endif
  }

  lua_config = config[config_filename]->globals()["config"];

  config_absorber = lua_config["absorber"].value_ptr<Absorber>();
  // Allow this to fail, we don't have aerosols if we happen to have a 
  // Rayleigh only atmosphere
  try {
    config_aerosol = lua_config["aerosol"].value_ptr<Aerosol>();
  } catch(const std::exception& e) {
    ;
  }
  config_atmosphere = lua_config["atmosphere"].value_ptr<RtAtmosphere>();
  config_state_vector = lua_config["state_vector"].value_ptr<StateVector>();
  config_pressure = lua_config["pressure"].value_ptr<Pressure>();
  config_instrument = lua_config["instrument"].value_ptr<Instrument>();
  config_spectral_window = lua_config["spec_win"].value_ptr<SpectralWindow>();
  config_initial_guess = 
    lua_config["initial_guess"].value_ptr<InitialGuess>();
  config_solver = lua_config["conn_solver"].value_ptr<ConnorSolver>();
  // Allow this to fail, we don't have ground if we happen to be looking
  // up.
  try {
    config_ground = lua_config["ground"].value_ptr<Ground>();
  } catch(const std::exception& e) {
    ;
  }
  config_temperature = lua_config["temperature"].value_ptr<Temperature>();
  config_spectrum_sampling = 
    lua_config["spec_samp"].value_ptr<SpectrumSampling>();
  config_error_analysis = 
    lua_config["error_analysis"].value_ptr<ErrorAnalysis>();
  config_level_1b = lua_config["l1b"].value_ptr<Level1bSampleCoefficient>();
  config_rt = lua_config["rt"].value_ptr<RadiativeTransfer>();
  config_forward_model = lua_config["forward_model"].value_ptr<ForwardModel>();
  config_observation = lua_config["observation"].value_ptr<Observation>();
  sv_initial.reference(config_initial_guess->initial_guess());
  config_state_vector->update_state(sv_initial);

}

void LuaConfigurationFixture::init_epsilon()
{
  epsilon.resize(config_state_vector->observer_claimed_size());
  epsilon = 1e-6;                  // Default
  epsilon(Range(0, 19)) = 1e-7;    // CO2 VMR
  epsilon(21) = 1e-3;              // Surface Pressure
  epsilon(22) = 1e-4;              // Temperature

  epsilon(Range(23, 34)) = 1e-8;   // Aerosol

  // For debugging epsilon values
  if (false) {
      for (int i=0 ; i < config_state_vector->observer_claimed_size(); i++) {
          std::cerr << i << ": " << config_state_vector->state_vector_name()(i) << " val = " << config_state_vector->state()(i) << ", epsilon = " << epsilon(i) << std::endl;
      }
  }
}

ConfigurationCoxmunkFixture::ConfigurationCoxmunkFixture(const std::string& Config_file) 
  : LuaConfigurationFixture(Config_file)
{
  // Shrink array from base constructor less elements but first many should be the same
  epsilon.resizeAndPreserve(110);
}


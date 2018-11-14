#ifndef CONFIGURATION_FIXTURE_H
#define CONFIGURATION_FIXTURE_H
#include "absorber.h"
#include "aerosol_optical.h"
#include "pressure_level_input.h"
#include "error_analysis.h"
#include "connor_solver.h"
#include "spectral_window.h"
#include "spectrum_sampling.h"
#include "level_1b_sample_coefficient.h"
#include "initial_guess.h"
#include "register_output_base.h"
#include "instrument.h"
#include "sample_grid.h"
#include "rt_atmosphere.h"
#include "radiative_transfer.h"
#include "lua_state.h"
#include "global_fixture.h"
#include "forward_model.h"
#include "observation.h"
#include <map>

namespace FullPhysics {
/****************************************************************//**
  There are a number of tests that need to use a 
  standard set of objects, similar to what we generate when running
  l2_fp. This is fairly expensive to create, so 
  this fixture creates one copy for all the tests and add handling for
  sharing it.

  Fixtures that are closely related to this one can derive from it and
  do things like change the config file to load (through the
  constructor argument), or individual values.
*******************************************************************/
class ConfigurationFixture: public GlobalFixture {
public:
  ConfigurationFixture(const std::string& Config_file = "config.lua");
  virtual ~ConfigurationFixture() 
  { config_state_vector->update_state(sv_initial); }

  boost::shared_ptr<Absorber> config_absorber;
  boost::shared_ptr<Aerosol> config_aerosol;
  boost::shared_ptr<RtAtmosphere> config_atmosphere;
  boost::shared_ptr<StateVector> config_state_vector;
  boost::shared_ptr<PressureLevelInput> config_pressure_level_input;
  boost::shared_ptr<ErrorAnalysis> config_error_analysis;
  boost::shared_ptr<ConnorSolver> config_solver;
  boost::shared_ptr<SpectralWindow> config_spectral_window;
  boost::shared_ptr<InitialGuess> config_initial_guess;
  std::vector<boost::shared_ptr<RegisterOutputBase> > config_register_output;
  boost::shared_ptr<Instrument> config_instrument;
  boost::shared_ptr<Pressure> config_pressure;
  boost::shared_ptr<Temperature> config_temperature;
  boost::shared_ptr<SpectrumSampling> config_spectrum_sampling;
  boost::shared_ptr<Level1bSampleCoefficient> config_level_1b;
  boost::shared_ptr<Ground> config_ground;
  boost::shared_ptr<ForwardModel> config_forward_model;
  boost::shared_ptr<Observation> config_observation;
  boost::shared_ptr<RadiativeTransfer> config_rt;

  /// Useful function that appears in a few tests, so collect here
  DoubleWithUnit ils_half_width(int Spec_index) const
  { return config_instrument->ils_half_width(Spec_index); }
  SpectralDomain lowres_grid(int Spec_index) const
  { return config_spectral_window->apply
      (config_instrument->pixel_spectral_domain(Spec_index), Spec_index); }
  SpectralDomain highres_grid(int Spec_index) const
  { return config_spectrum_sampling->spectral_domain
      (Spec_index, lowres_grid(Spec_index), ils_half_width(Spec_index)); }
    
  /// This is an epsilon that can be used to generate finite
  /// difference Jacobians for unit tests.
  blitz::Array<double, 1> epsilon;

  LuabindObject lua_config;
private:
  static std::map<std::string, boost::shared_ptr<LuaState> > config;
  blitz::Array<double,1> sv_initial;
};

/****************************************************************//**
  There are a number of tests that need to use a 
  standard set of objects, similar to what we generate when running
  l2_fp. This is fairly expensive to create, so 
  this fixture creates one copy for all the tests and add handling for
  sharing it. This version reads config_coxmunk.lua
*******************************************************************/
class ConfigurationCoxmunkFixture: public ConfigurationFixture {
public:
  ConfigurationCoxmunkFixture(const std::string& Config_file = "config_coxmunk.lua");
  virtual ~ConfigurationCoxmunkFixture() {}
};

/****************************************************************//**
  There are a number of tests that need to use a 
  standard set of objects, similar to what we generate when running
  l2_fp. This is fairly expensive to create, so 
  this fixture creates one copy for all the tests and add handling for
  sharing it. This version reads config_coxmunk+lamb.lua
*******************************************************************/
class ConfigurationCoxmunkPlusLambertianFixture: public ConfigurationFixture {
public:
  ConfigurationCoxmunkPlusLambertianFixture()
    : ConfigurationFixture("config_coxmunk+lamb.lua") {}
  virtual ~ConfigurationCoxmunkPlusLambertianFixture() {}
};

/****************************************************************//**
  There are a number of tests that need to use a 
  standard set of objects, similar to what we generate when running
  l2_fp. This is fairly expensive to create, so 
  this fixture creates one copy for all the tests and add handling for
  sharing it. This version reads config_brdf_veg.lua
*******************************************************************/
class ConfigurationBrdfVegFixture: public ConfigurationFixture {
public:
  ConfigurationBrdfVegFixture()
    : ConfigurationFixture("config_brdf_veg.lua") {}
  virtual ~ConfigurationBrdfVegFixture() {}
};

/****************************************************************//**
  There are a number of tests that need to use a 
  standard set of objects, similar to what we generate when running
  l2_fp. This is fairly expensive to create, so 
  this fixture creates one copy for all the tests and add handling for
  sharing it. This version reads config_brdf_soil.lua
*******************************************************************/
class ConfigurationBrdfSoilFixture: public ConfigurationFixture {
public:
  ConfigurationBrdfSoilFixture()
    : ConfigurationFixture("config_brdf_soil.lua") {}
  virtual ~ConfigurationBrdfSoilFixture() {}
};

}
#endif

#ifndef LUA_CONFIGURATION_FIXTURE_H
#define LUA_CONFIGURATION_FIXTURE_H

#include "configuration_fixture.h"
#include "error_analysis.h"
#include <map>

namespace FullPhysics {
/****************************************************************//**
  This fixture loads ConfigurationFixture objects by loading
  a Lua configuration.

  This ConfigurationFixture is in the process of being deprecated
  in preference for serialized objects.
*******************************************************************/
class LuaConfigurationFixture: public ConfigurationFixture {
public:

  LuaConfigurationFixture(const std::string& Config_file = "config.lua");

  virtual ~LuaConfigurationFixture() 
  { config_state_vector->update_state(sv_initial); }

  boost::shared_ptr<ErrorAnalysis> config_error_analysis;

  LuabindObject lua_config;
protected:
  virtual void init_variables();
  virtual void init_epsilon();
private:
  std::string config_filename;
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
class ConfigurationCoxmunkFixture: public LuaConfigurationFixture {
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
class ConfigurationCoxmunkPlusLambertianFixture: public LuaConfigurationFixture {
public:
  ConfigurationCoxmunkPlusLambertianFixture()
    : LuaConfigurationFixture("config_coxmunk+lamb.lua") {}
  virtual ~ConfigurationCoxmunkPlusLambertianFixture() {}
};

/****************************************************************//**
  There are a number of tests that need to use a 
  standard set of objects, similar to what we generate when running
  l2_fp. This is fairly expensive to create, so 
  this fixture creates one copy for all the tests and add handling for
  sharing it. This version reads config_brdf_veg.lua
*******************************************************************/
class ConfigurationBrdfVegFixture: public LuaConfigurationFixture {
public:
  ConfigurationBrdfVegFixture()
    : LuaConfigurationFixture("config_brdf_veg.lua") {}
  virtual ~ConfigurationBrdfVegFixture() {}
};

/****************************************************************//**
  There are a number of tests that need to use a 
  standard set of objects, similar to what we generate when running
  l2_fp. This is fairly expensive to create, so 
  this fixture creates one copy for all the tests and add handling for
  sharing it. This version reads config_brdf_soil.lua
*******************************************************************/
class ConfigurationBrdfSoilFixture: public LuaConfigurationFixture {
public:
  ConfigurationBrdfSoilFixture()
    : LuaConfigurationFixture("config_brdf_soil.lua") {}
  virtual ~ConfigurationBrdfSoilFixture() {}
};

/****************************************************************//**
  There are a number of tests that need to use a 
  standard set of objects, similar to what we generate when running
  l2_fp. This is fairly expensive to create, so 
  this fixture creates one copy for all the tests and add handling for
  sharing it. This version reads config_two_broadener.lua
*******************************************************************/
class ConfigurationTwoBroadener: public LuaConfigurationFixture {
public:
  ConfigurationTwoBroadener()
    : LuaConfigurationFixture("config_two_broadener.lua") {}
  virtual ~ConfigurationTwoBroadener() {}
};
}
#endif

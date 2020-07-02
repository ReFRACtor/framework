# We are moving away from the Lua configuration, but for now there isn't
# any easy way to generate a lot of the more complicated test object (e.g.
# a forward model). We'll try to restrict this stuff to the fixtures here,
# so there is one place to update when we can move away from Lua.

from refractor import framework as rf
from refractor import read_shelve, write_shelve
import os
import sys
import pytest

# Have unit test LUA_PATH set.

os.environ["LUA_PATH"] = os.path.abspath(os.path.dirname(__file__) + "/../../../test/unit/data/lua") + "/?.lua;" + os.path.abspath(os.path.dirname(__file__) + "/../../../input/common/config") + "/?.lua"

serialized_test_data = os.path.abspath(os.path.dirname(__file__)) + "/"

@pytest.fixture(scope="function")
def config_ls():
    return rf.LuaState.load_file( os.path.abspath(os.path.dirname(__file__) + "/../../../test/unit/data/lua") + "/config.lua")


@pytest.fixture(scope="function")
def config_forward_model(config_ls):
    '''A forward model example. Currently comes from Lua, but we may
    change the source in the future'''
    return config_ls.globals.config.forward_model

@pytest.fixture(scope="function")
def sample_absorber(config_ls):
    '''Sample Absorber, to use in tests that need one. 

    We save the data in a serialization file. The file gets regenerated
    if it is missing, so if you need to regenerate just delete the file
    and run the tests.
    
    This currently comes from the old Lua code, just for convenience. We
    probably want to come up with a different method (e.g, directly create
    this), since we are trying to phase out the Lua. But for now, use it 
    because it is available.'''
    return config_ls.globals.config.absorber

@pytest.fixture(scope="function")
def sample_temperature(config_ls):
    '''Sample temperature, to use in tests that need one. 

    This currently comes from the old Lua code, just for convenience. We
    probably want to come up with a different method (e.g, directly create
    this), since we are trying to phase out the Lua. But for now, use it 
    because it is available.'''
    return config_ls.globals.config.temperature

@pytest.fixture(scope="function")
def sample_pressure(config_ls):
    '''Sample pressure, to use in tests that need one. 

    This currently comes from the old Lua code, just for convenience. We
    probably want to come up with a different method (e.g, directly create
    this), since we are trying to phase out the Lua. But for now, use it 
    because it is available.'''
    return config_ls.globals.config.pressure

@pytest.fixture(scope="function")
def sample_aerosol(config_ls):
    '''Sample aerosol, to use in tests that need one. 

    This currently comes from the old Lua code, just for convenience. We
    probably want to come up with a different method (e.g, directly create
    this), since we are trying to phase out the Lua. But for now, use it 
    because it is available.'''
    return config_ls.globals.config.aerosol

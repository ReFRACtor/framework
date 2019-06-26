# We are moving away from the Lua configuration, but for now there isn't
# any easy way to generate a lot of the more complicated test object (e.g.
# a forward model). We'll try to restrict this stuff to the fixtures here,
# so there is one place to update when we can move away from Lua.

from refractor import framework as rf
import os
import sys
import pytest

# Have unit test LUA_PATH set.

os.environ["LUA_PATH"] = os.path.abspath(os.path.dirname(__file__) + "/../../../test/unit/data/lua") + "/?.lua;" + os.path.abspath(os.path.dirname(__file__) + "/../../../input/common/config") + "/?.lua"

@pytest.yield_fixture(scope="function")
def config_ls():
    return rf.LuaState.load_file( os.path.abspath(os.path.dirname(__file__) + "/../../../test/unit/data/lua") + "/config.lua")


@pytest.yield_fixture(scope="function")
def config_forward_model(config_ls):
    '''A forward model example. Currently comes from Lua, but we may
    change the source in the future'''
    return config_ls.globals.config.forward_model

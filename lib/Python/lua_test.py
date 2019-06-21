from .test_support import *
import os

# A callback that takes no arguments
class MyCallback0(rf.LuaCallback):
    def __init__(self, ls):
        rf.LuaCallback.__init__(self, ls)
        self.ls = ls
        
    def call(self, obj1,obj2,obj3,obj4,obj5,obj6,obj7,obj8,obj9,obj10):
        return rf.LuabindObject(self.ls, 1)
        
# A callback that takes one arguments
class MyCallback1(rf.LuaCallback):
    def __init__(self, ls):
        rf.LuaCallback.__init__(self, ls)
        self.ls = ls
                
    def call(self, obj1,obj2,obj3,obj4,obj5,obj6,obj7,obj8,obj9,obj10):
        obj1.v3 = "hi there"
        return rf.LuabindObject(self.ls, 2)

# HeritageFile doesn't exist any longer, so skip test    
@skip    
def test_read_and_write(lua_state):
    '''Test basic reading and writing of object to Lua'''
    lua_state.run("test_var = {}\n")
    test_var = lua_state.globals.test_var
    test_var.v1 = "hi there"
    test_var.v2 = 1
    test_var.v3 = 2.0
    t = HeritageFile(unit_test_data + "heritage_file_test.run")
    assert t.value_int("ALGORITHMS/points_sun") == 10000
    test_var.test_obj = t
    assert lua_state.globals.test_var.v1 == "hi there"
    assert lua_state.globals.test_var.v2 == 1
    assert lua_state.globals.test_var.v3 == 2.0
    assert lua_state.globals.test_var.test_obj.value_int("ALGORITHMS/points_sun") == 10000

    test_list = ["a", "b", "abc"]
    test_var.test_list = test_list
    for l_idx, l_val in enumerate(test_list):
        assert lua_state.globals.test_var.test_list[l_idx+1] == test_list[l_idx]

    test_dict = {"k1": "abc", "k2": "cde", 1: "blah"}
    test_var.test_dict = test_dict
    for l_key, l_val in list(test_dict.items()):
        assert lua_state.globals.test_var.test_dict[l_key] == test_dict[l_key]

@skip
def test_config(lua_state):
    '''Test reading a config file'''
    lua_state.do_file(unit_test_data + "/lua/config.lua")
    assert lua_state.globals.config.atmosphere.pressure.surface_pressure.value.value == 96716.6249

def test_callback_object(lua_state):
    '''Test a callback function object. Note you don't normally use this 
    directly, but rather use the callback tested in the next function'''
    lua_state.run("test_var = {}\n")
    test_var = lua_state.globals.test_var
    test_var.f0 = MyCallback0(lua_state)
    test_var.f1 = MyCallback1(lua_state)
    global test_v
    test_v = 0
    assert test_v == 0
    lua_state.run("test_v = test_var.f0()")
    assert lua_state.globals.test_v == 1
    lua_state.run("test_v = test_var:f1()")
    assert lua_state.globals.test_v == 2
    assert test_var.v3 == "hi there"

def f0():
    global test_v
    test_v = 1

def f1(v):
    global test_v
    test_v = v

def f2(v1, v2):
    global test_v
    test_v = v1 + v2

def f3(v1, v2, v3):
    global test_v
    test_v = v1 + v2 + v3

def f3_return(v1, v2, v3):
    return v1 + v2 + v3

def test_callback(lua_state):
    '''Test passing any callable object to Lua, to make sure the callbacks 
    work'''
    g = lua_state.globals
    g.f0 = f0
    g.f1 = f1
    g.f2 = f2
    g.f3 = f3
    g.f3_return = f3_return
    global test_v
    test_v = 0
    assert test_v == 0
    lua_state.run("f0(nil,nil,nil,nil,nil,nil,nil,nil,nil)")
    assert test_v == 1
    lua_state.run("f1(3.5)")
    assert test_v == 3.5
    lua_state.run("f2(3, 4)")
    assert test_v == 3 + 4
    lua_state.run("f3(3, 4, 5)")
    assert test_v == 3 + 4 + 5
    lua_state.run("val = f3_return(3, 4, 5)")
    assert g.val == 3 + 4 + 5

def test_luafunc(lua_state):
    '''Test calling a Lua function.'''
    lua_state.run('''
                    function test_func(v1, v2)
                       return v1 + v2 * 2
                    end''')
    tf = lua_state.globals.test_func
    assert tf(1, 2) == 1 + 2 * 2
    lua_state.run('''
                    function test_func2()
                       return "blah"
                    end''')
    tf2 = lua_state.globals.test_func2
    assert tf2() == "blah"
        

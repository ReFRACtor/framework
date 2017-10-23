import numpy as np
from functools import partial

import refractor.factory.creator as creator
import refractor.factory.param as param
from refractor.factory import process_config

from refractor import framework as rf

def ParamReturnCreator(param_type, **kwargs):
    class ReturnCreatorHelper(creator.base.Creator):
        val = param_type(**kwargs)
        def create(self, **kwargs):
            return self.param("val")

    return ReturnCreatorHelper

class AddCreator(creator.base.Creator):
    
    x = param.Scalar()
    y = param.Scalar()
    
    def create(self, **kwargs):
        x = self.param("x")
        y = self.param("y")
        return x + y

def test_scalar_param():
    
    # Check that scalar params are correctly checked
    config_def = {
        'creator': ParamReturnCreator(param.Scalar),
        'val': 10,
    }
 
    config_inst = process_config(config_def)

    assert config_inst == 10

    # Check that non scalars cause errors
    config_def['val'] = np.arange(1,6)

    try:
        config_inst = process_config(config_def)
    except param.ParamError:
        assert True
    else:
        assert False

def test_scalar_with_type():

    # Check that scalar params are correctly checked
    config_def = {
        'creator': ParamReturnCreator(partial(param.Scalar, dtype=str)),
        'val': "10",
    }
 
    config_inst = process_config(config_def)

    assert config_inst == "10"

    # Check that non scalars cause errors
    config_def['val'] = 10

    try:
        config_inst = process_config(config_def)
    except param.ParamError:
        assert True
    else:
        assert False

def test_array_param():
    
    # Check that arrays are correctly handled
    config_def = {
        'creator': ParamReturnCreator(param.Array),
        'val': np.arange(1,6),
    }
 
    config_inst = process_config(config_def)

    assert np.all(config_inst == np.array([1,2,3,4,5]))

    # Check that non arrays cause errors
    config_def['val'] = 10

    try:
        config_inst = process_config(config_def)
    except param.ParamError:
        assert True
    else:
        assert False

    config_def['val'] = [1,2,3,4,5]

    try:
        config_inst = process_config(config_def)
    except param.ParamError:
        assert True
    else:
        assert False

def test_array_with_unit_param():

    # Check that arrays are correctly handled
    config_def = {
        'creator': ParamReturnCreator(param.ArrayWithUnit),
        'val': rf.ArrayWithUnit_double_1(np.arange(1,6), "m"),
    }

    config_inst = process_config(config_def)

 
def test_iterable_param():

    # Check that iterables are correctly handled
    config_def = {
        'item1': {
            'creator': ParamReturnCreator(param.Iterable),
            'val': [1,2,3,4,5]
        },
        'item2': {
            'creator': ParamReturnCreator(param.Iterable, val_type=np.int),
            'val': [1,2,3,4,5]
        },
    }
 
    config_inst = process_config(config_def)

    assert np.all(np.array(config_inst['item1']) == np.arange(1,6))
    assert np.all(np.array(config_inst['item2']) == np.arange(1,6))

    # Check arrays are considered iterable
    config_def['item1']['val'] = np.arange(1,6)

    config_inst = process_config(config_def)

    assert np.all(np.array(config_inst['item1']) == np.arange(1,6))

def test_instanceof_param():

    class TestClass(object):
        pass

    # Check that iterables are correctly handled
    config_def = {
        'creator': ParamReturnCreator(param.InstanceOf, val_type=TestClass),
        'val': TestClass(),
    }
 
    config_inst = process_config(config_def)

    assert isinstance(config_inst, TestClass)

    # Force error on type
    config_def['val'] = 5

    try:
        config_inst = process_config(config_def)
    except param.ParamError:
        assert True
    else:
        assert False

def test_callable():

    config_def = {
        'creator': ParamReturnCreator(param.Scalar),
        'val': lambda **kwargs: 10,
    }

    config_inst = process_config(config_def)

    assert config_inst == 10

def test_nested():

    config_def = {
        'item1': {
            'creator': ParamReturnCreator(param.Scalar),
            'val': 5,
        },
        'item2': {
            'creator': ParamReturnCreator(param.Scalar),
            'val': 10,
        },
        'item3': {
            'creator': ParamReturnCreator(param.Scalar),
            'val': {
                'creator': ParamReturnCreator(param.Scalar),
                'val': 15,
            },
        },
    }

    config_inst = process_config(config_def)

    assert config_inst['item1'] == 5
    assert config_inst['item2'] == 10
    assert config_inst['item3'] == 15

def test_common_store():

    config_def = {
        'order': ['common', 'use_common'],
        'common': {
            'creator': creator.base.SaveToCommon,
            'x': 5,
            'y': 6,
        },
        'use_common': {
            'creator': AddCreator,
        },
    }

    config_inst = process_config(config_def)

    assert config_inst['use_common'] == 11

def test_param_choice():

    class ChoiceCreator(creator.base.Creator):
        some_val = param.Choice(param.Scalar(int), param.Scalar(float))

        def create(self, **kwargs):
            return self.param("some_val")

    config_def = { 
        'item1': { 
            'creator': ChoiceCreator,
            'some_val': 10,
        },
        'item2': { 
            'creator': ChoiceCreator,
            'some_val': 5.0,
        },
    }

    config_inst = process_config(config_def)

    assert config_inst['item1'] == 10
    assert config_inst['item2'] == 5.0

    config_def['item1']['some_val'] = "ABC"

    try:
        config_inst = process_config(config_def)
    except param.ParamError:
        assert True
    else:
        assert False

def test_bound_params():

    class BoundParamCreator(creator.base.Creator):
        some_val = param.Choice(param.Scalar(int), param.Scalar(int))

        def create(self, **kwargs):
            return self.some_val()

    config_def = { 
        'item1': { 
            'creator': BoundParamCreator,
            'some_val': 10,
        },
        'item2': { 
            'creator': BoundParamCreator,
            'some_val': 5,
        },
    }

    config_inst = process_config(config_def)

    assert config_inst['item1'] == 10
    assert config_inst['item2'] == 5

def test_object_vector():

    config_def = {
        'item1': {
            'creator': ParamReturnCreator(param.ObjectVector),
            'val': rf.vector_double(),
        },
        'item2': {
            'creator': ParamReturnCreator(param.ObjectVector, vec_type="double"),
            'val': rf.vector_double(),
        }
    }

    config_inst = process_config(config_def)

    assert isinstance(config_inst['item1'], rf.vector_double)
    assert isinstance(config_inst['item2'], rf.vector_double)

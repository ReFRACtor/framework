import numpy as np

import full_physics.factory.creator as creator
import full_physics.factory.param as param
from full_physics.factory import process_config

def ParamReturnCreator(param_type, **kwargs):
    class ReturnCreatorHelper(creator.base.Creator):
        val = param_type(**kwargs)
        def create(self):
            return self.param("val")

    return ReturnCreatorHelper

class AddCreator(creator.base.Creator):
    
    x = param.Scalar()
    y = param.Scalar()
    
    def create(self):
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

def test_iterable_param():

    # Check that iterables are correctly handled
    config_def = {
        'creator': ParamReturnCreator(param.Iterable),
        'val': [1,2,3,4,5]
    }
 
    config_inst = process_config(config_def)

    assert np.all(np.array(config_inst) == np.arange(1,6))

    # Check arrays are considered iterable
    config_def['val'] = np.arange(1,6)

    config_inst = process_config(config_def)

    assert np.all(np.array(config_inst) == np.arange(1,6))

def test_instanceof_param():

    class TestClass(object):
        pass

    # Check that iterables are correctly handled
    config_def = {
        'creator': ParamReturnCreator(param.InstanceOf, cls_type=TestClass),
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
        'val': lambda creator: 10,
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

# Provides fixtures with sample data that can be utilized by test cases

import refractor.framework as rf
import os
import pytest

serialized_test_data_dir = os.path.abspath(os.path.dirname(__file__)) + "/../../../test/unit/data/in/configuration_fixture/"

@pytest.fixture(scope="function")
def lambertian_objects():
    objs = rf.read_shelve(os.path.join(serialized_test_data_dir, "lambertian_example_config.bin.gz"))
    return objs

@pytest.fixture(scope="function")
def sample_forward_model(lambertian_objects):
    '''A forward model example.'''

    return lambertian_objects.forward_model

@pytest.fixture(scope="function")
def sample_solar_model(sample_forward_model):
    '''A solar model sample example object'''

    return sample_forward_model.find_subobject_of_type(rf.SolarModel)

@pytest.fixture(scope="function")
def sample_absorber(sample_forward_model):
    '''Sample Absorber, to use in tests that need one.'''

    return sample_forward_model.find_subobject_of_type(rf.Absorber)

@pytest.fixture(scope="function")
def sample_temperature(sample_forward_model):
    '''Sample temperature, to use in tests that need one.'''

    return sample_forward_model.find_subobject_of_type(rf.Temperature)

@pytest.fixture(scope="function")
def sample_atmosphere(sample_forward_model):
    '''Sample atmosphere, to use in tests that need one.'''

    return sample_forward_model.find_subobject_of_type(rf.RtAtmosphere)

@pytest.fixture(scope="function")
def sample_pressure(sample_forward_model):
    '''Sample pressure, to use in tests that need one.'''

    return sample_forward_model.find_subobject_of_type(rf.Pressure)

@pytest.fixture(scope="function")
def sample_aerosol(sample_forward_model):
    '''Sample aerosol, to use in tests that need one.'''

    return sample_forward_model.find_subobject_of_type(rf.Aerosol)

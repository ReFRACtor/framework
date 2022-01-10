#!/usr/bin/env python3

import os
import sys
import logging

from refractor import framework as rf
from refractor.framework import load_config_module, find_config_function
from refractor.framework.factory import process_config
from refractor.framework import write_shelve, GenericObjectMap

logging.basicConfig(level=logging.DEBUG)

config_dir = os.path.realpath(os.path.join(os.path.dirname(__file__), "../config"))
input_dir = os.path.realpath(os.path.join(os.path.dirname(__file__), "../data/in"))

config_fn = os.path.join(config_dir, "example_base_config.py")
serialized_fn = os.path.join(input_dir, "configuration_fixture/example_base_config.bin.gz")

logging.info(f"Loading configuration: {config_fn}")

config_module = load_config_module(config_fn)
config_func = find_config_function(config_module)

config_def = config_func()

config_inst = process_config(config_def)

# Update state vector with initial guess so that individual SubStateVectorObserver objects
# have a value for their sv_full values and the state_vector state itself is non zero
config_inst.retrieval.state_vector.update_state(config_inst.retrieval.initial_guess)

logging.info(f"Writing serialized representation: {serialized_fn}")

# Copy objects from Python config into a serializable object in the form
# used by the ConfigurationFixture classes
obj_map = GenericObjectMap()

# Copy base objects
config_object_names = [
    "atmosphere",
    "instrument",
    "spectrum_sampling",
    "forward_model",
]

for config_name in config_object_names:
    obj_map[config_name] = config_inst[config_name]

# Copy atmosphere objects
config_object_names = [
    "absorber",
    "aerosol",
    "pressure",
    "temperature",
    "ground",
]

for config_name in config_object_names:
    obj_map[config_name] = getattr(config_inst.atmosphere, config_name)

# Copy retrieval objects
config_object_names = [
    "state_vector",
    "solver",
]

for config_name in config_object_names:
    obj_map[config_name] = config_inst.retrieval[config_name]

# Initial guess in Python configuration is reperesented by an array, wrap in an InitialGuess object
# to meet ConfigurationFixture interface
initial_guess_builder = rf.InitialGuessValue()
initial_guess_builder.apriori = config_inst.retrieval.initial_guess
initial_guess_val = rf.CompositeInitialGuess()
initial_guess_val.add_builder(initial_guess_builder)
obj_map["initial_guess"] = initial_guess_val
    

# Remaining objects
obj_map["spectral_window"] = config_inst.spec_win
obj_map["level_1b"] = config_inst.input.l1b
obj_map["observation"] = config_inst.retrieval.solver.cost_function.observation
obj_map["rt"] = config_inst.radiative_transfer

write_shelve(serialized_fn, obj_map)

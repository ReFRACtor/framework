#!/usr/bin/env python3

import os
import sys
import logging

from refractor import framework as rf
from refractor.framework import load_config_module, find_config_function
from refractor.framework.factory import creator, process_config
from refractor.framework import write_shelve, GenericObjectMap

logging.basicConfig(level=logging.INFO)

config_dir = os.path.realpath(os.path.join(os.path.dirname(__file__), "../config"))
input_dir = os.path.realpath(os.path.join(os.path.dirname(__file__), "../data/in"))

# Define this environment variable so that it can be used within configuration to set paths that get expended out during unit testing execution
os.environ['abs_top_srcdir'] = os.path.realpath(os.path.join(os.path.dirname(__file__), "../../.."))

class SerializeConfigObjects(object):

    def __init__(self, config_fn):

        self.config_fn = config_fn

        self._load_config_definition()
        self._instantiate_config()
        self._extract_objects()

    def _load_config_definition(self):
        """Load the configuration definition structure"""

        logging.info(f"Loading configuration: {self.config_fn}")

        config_module = load_config_module(self.config_fn)
        config_func = find_config_function(config_module)

        self.config_def = config_func()

    def _instantiate_config(self):
        """Loads config definition and processes it into an instantiated config"""

        logging.info(f"Processing configuration")

        self.config_inst = process_config(self.config_def)

        # Update state vector with initial guess so that individual SubStateVectorObserver objects
        # have a value for their sv_full values and the state_vector state itself is non zero
        self.config_inst.retrieval.state_vector.update_state(self.config_inst.retrieval.initial_guess)

    def _extract_objects(self):
        """Copy objects from Python config into a serializable object in the form
           used by the ConfigurationFixture classes"""

        logging.info("Extracting objects to serialize")

        self.obj_map = GenericObjectMap()

        # Copy base objects
        config_object_names = [
            "atmosphere",
            "instrument",
            "spectrum_sampling",
            "forward_model",
        ]

        for config_name in config_object_names:
            self.obj_map[config_name] = self.config_inst[config_name]

        # Copy atmosphere objects
        config_object_names = [
            "absorber",
            "aerosol",
            "pressure",
            "temperature",
            "ground",
        ]

        for config_name in config_object_names:
            self.obj_map[config_name] = getattr(self.config_inst.atmosphere, config_name)

        # Copy retrieval objects
        config_object_names = [
            "state_vector",
            "solver",
        ]

        for config_name in config_object_names:
            self.obj_map[config_name] = self.config_inst.retrieval[config_name]

        # Initial guess in Python configuration is reperesented by an array, wrap in an InitialGuess object
        # to meet ConfigurationFixture interface
        initial_guess_builder = rf.InitialGuessValue()
        initial_guess_builder.apriori = self.config_inst.retrieval.initial_guess
        initial_guess_builder.apriori_covariance = self.config_inst.retrieval.covariance
        initial_guess_val = rf.CompositeInitialGuess()
        initial_guess_val.add_builder(initial_guess_builder)
        self.obj_map["initial_guess"] = initial_guess_val
            

        # Remaining objects
        self.obj_map["spectral_window"] = self.config_inst.spec_win
        self.obj_map["level_1b"] = self.config_inst.input.l1b
        self.obj_map["observation"] = self.config_inst.retrieval.solver.cost_function.observation
        self.obj_map["rt"] = self.config_inst.radiative_transfer

    def dump(self, serialized_fn):
        """Serialized objects into an output file"""

        logging.info(f"Writing serialized objects: {serialized_fn}")
        write_shelve(serialized_fn, self.obj_map)

class SerializeWithFluorescence(SerializeConfigObjects):

    def _load_config_definition(self):
        super()._load_config_definition()
        
        self.config_def["fluorescence"] = creator.util.ObjectCapture(rf.FluorescenceEffect)
        self.config_def["order"].append("fluorescence")

    def _extract_objects(self):
        super()._extract_objects()

        self.obj_map["fluorescence"] = self.config_inst.fluorescence
 
def serialize_config(config_fn, serialized_fn):

    ser = SerializeConfigObjects(config_fn)
    ser.dump(serialized_fn)

def serialize_fluor_config(config_fn, serialized_fn):

    ser = SerializeWithFluorescence(config_fn)
    ser.dump(serialized_fn)

serialize_config(os.path.join(config_dir, "lambertian_example_config.py"),
                 os.path.join(input_dir, "configuration_fixture/lambertian_example_config.bin.gz"))

serialize_config(os.path.join(config_dir, "coxmunk_example_config.py"),
                 os.path.join(input_dir, "configuration_fixture/coxmunk_example_config.bin.gz"))

serialize_config(os.path.join(config_dir, "coxmunk_lambertian_example_config.py"),
                 os.path.join(input_dir, "configuration_fixture/coxmunk_lambertian_example_config.bin.gz"))

serialize_config(os.path.join(config_dir, "brdf_veg_example_config.py"),
                 os.path.join(input_dir, "configuration_fixture/brdf_veg_example_config.bin.gz"))

serialize_config(os.path.join(config_dir, "brdf_soil_example_config.py"),
                 os.path.join(input_dir, "configuration_fixture/brdf_soil_example_config.bin.gz"))

serialize_config(os.path.join(config_dir, "two_broadener_example_config.py"),
                 os.path.join(input_dir, "configuration_fixture/two_broadener_example_config.bin.gz"))

serialize_fluor_config(os.path.join(config_dir, "fluorescence_example_config.py"),
                       os.path.join(input_dir, "configuration_fixture/fluorescence_example_config.bin.gz"))

serialize_fluor_config(os.path.join(config_dir, "rayleigh_only_example_config.py"),
                       os.path.join(input_dir, "configuration_fixture/rayleigh_only_example_config.bin.gz"))

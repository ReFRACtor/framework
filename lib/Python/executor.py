import os
import sys
import logging
import pprint

from netCDF4 import Dataset

from .config import load_config_module, find_config_function, find_strategy_function
from .factory import process_config, creator
from . import framework as rf

from .output.radiance import ForwardModelRadianceOutput, ObservationRadianceOutput
from .output.solver import SolverIterationOutput
from .output.state_vector import StateVectorOutputRetrieval, StateVectorOutputSimulation

logger = logging.getLogger(__name__)

def ObjectCapture(capture_class):
    "Generate a Creator that watches for an object to be emitted elsewhere then stores it internally to use as its return object"

    class ObjectCaptureCreator(creator.base.Creator):
        def __init__(self, *vargs, **kwargs):
            super().__init__(*vargs, **kwargs)

            self.register_to_receive(capture_class)
            self.captured_object = None

        def receive(self, rec_obj):
            self.captured_object = rec_obj

        def create(self, **kwargs): 
            return self.captured_object

    return ObjectCaptureCreator

class StrategyExecutor(object):
    
    def __init__(self, config_filename, output_filename=None, strategy_filename=None, strategy_list=None):

        logger.debug("Loading configuration from {}".format(config_filename))
        self.config_module = load_config_module(config_filename)

        if strategy_filename is not None:
            logger.debug("Loading strategy from {}".format(strategy_filename))
            strategy_module = load_config_module(strategy_filename)
            self.strategy_list = find_strategy_function(strategy_module)()
        elif strategy_list is not None:
            self.strategy_list = strategy_list
        else:
            self.strategy_list = [ {} ]

        logger.debug("{} strategy steps loaded".format(len(self.strategy_list)))

        if output_filename is not None:
            self.output = Dataset(output_filename, "w")
        else:
            self.output = None

    def config_definition(self, **strategy_keywords):

        config_func = find_config_function(self.config_module)

        # Augment loaded configuration to help capture objects required for ouput
        config_def = {
            'creator': creator.base.ParamPassThru, 
            'order': ['file_config', 'forward_model', 'solver', 'state_vector', 'l1b'],
            'file_config': config_func(**strategy_keywords),

            # Capture certain objects without depending on the structure of the configuration
            'forward_model': { 'creator': ObjectCapture(rf.ForwardModel) },
            'solver': { 'creator': ObjectCapture(rf.IterativeSolver) },
            'state_vector': { 'creator': ObjectCapture(rf.StateVector) },
            'l1b': { 'creator': ObjectCapture(rf.Level1b) },
        }

        return config_def

    def config_instance(self, **strategy_keywords):

        logger.debug("Loading configration")

        config_def = self.config_definition(**strategy_keywords)
        config_inst = process_config(config_def)

        logger.debug("Configuration processing complete")

        return config_inst

    def attach_logging(self, config_inst):

        iter_log = rf.SolverIterationLog(config_inst['state_vector'])
        config_inst['solver'].add_observer_and_keep_reference(iter_log)

    def attach_output(self, config_inst, step_index):

        if self.output is None:
            return

        rad_out = ForwardModelRadianceOutput(self.output, step_index, config_inst['solver'])
        config_inst['forward_model'].add_observer_and_keep_reference(rad_out)

        if config_inst['l1b'] is not None and config_inst['solver'] is not None:
            obs_out = ObservationRadianceOutput(self.output, step_index, config_inst['l1b'], config_inst['forward_model'])
            config_inst['solver'].add_observer_and_keep_reference(obs_out)

        if config_inst['state_vector'] is not None and config_inst['solver'] is not None:
            solver_out = SolverIterationOutput(self.output, step_index)
            config_inst['solver'].add_observer_and_keep_reference(solver_out)

            sv_out = StateVectorOutputRetrieval(self.output, step_index, config_inst['state_vector'])
            config_inst['solver'].add_observer_and_keep_reference(sv_out)

        elif config_inst['state_vector'] is not None and config_inst['solver'] is None:
            sv_out = StateVectorOutputSimulation(self.output, step_index, config_inst['state_vector'])
 
    def run_solver(self, config_inst, step_index=None):

        self.attach_logging(config_inst)
        self.attach_output(config_inst, step_index)

        if config_inst['solver'] is None:
            raise Exception("Solver object was not defined in configuration")

        if config_inst['state_vector'] is None:
            raise Exception("StateVector object was not defined in configuration")
 
        config_inst['solver'].solve()
        config_inst['state_vector'].clear_observers()

    def run_forward_model(self, config_inst, step_index=None):

        file_inst = config_inst['file_config']
        if config_inst['state_vector'] is not None and 'retrieval' in file_inst and 'initial_guess' in file_inst['retrieval']:
            # Necessary for jacobians to work. Normally this step is done in the solver creator
            # Do this before attaching output so that state vector can be output correctly
            logger.debug("Enabling jacobians for forward model simulation")
            config_inst['state_vector'].update_state(file_inst['retrieval']['initial_guess'])

        self.attach_output(config_inst, step_index)

        config_inst['forward_model'].radiance_all()

    def execute_retrieval(self, step_indexes=None):

        if step_indexes is None or len(step_indexes) == 0:
            step_indexes = range(len(self.strategy_list))
        else:
            logger.debug("Using {} steps".format(len(step_indexes)))

        for step_index in step_indexes:
            if step_index < 0 or step_index >= len(self.strategy_list):
                raise Exception("Invalid step index: {}".format(step_index))

            logger.info("Executing retrieval step #{}".format(step_index+1))

            step_keywords = self.strategy_list[step_index]

            logger.debug("Strategy step options:")
            logger.debug(pprint.pformat(step_keywords, indent=2))

            config_inst = self.config_instance(**step_keywords)

            self.run_solver(config_inst, step_index)

            del config_inst

        logger.info("Retrieval execution finished")

    def execute_simulation(self, step_indexes=None):

        if step_indexes is None or len(step_indexes) == 0:
            step_indexes = range(len(self.strategy_list))
        else:
            logger.debug("Using {} steps".format(len(step_indexes)))

        for step_index in step_indexes:
            if step_index < 0 or step_index >= len(self.strategy_list):
                raise Exception("Invalid step index: {}".format(step_index))

            logger.info("Executing simulation for step #{}".format(step_index+1))

            step_keywords = self.strategy_list[step_index]
            
            config_inst = self.config_instance(**step_keywords)

            self.run_forward_model(config_inst, step_index)

            del config_inst

        logger.info("Simulation execution finished")

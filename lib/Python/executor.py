import os
import sys
import logging
import pprint

import numpy as np
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

        # Location where covariances data should be saved
        self.covariance_storage = {}

    def config_definition(self, **strategy_keywords):

        config_func = find_config_function(self.config_module)


        # Augment loaded configuration to help capture objects required for ouput
        config_def = {
            'creator': creator.base.ParamPassThru, 
            'order': ['file_config'],
            'file_config': config_func(**strategy_keywords),
        }

        # Capture certain objects without depending on the structure of the configuration
        captured_objects = {
            'forward_model': ObjectCapture(rf.ForwardModel),
            'solver': ObjectCapture(rf.IterativeSolver),
            'state_vector': ObjectCapture(rf.StateVector),
            'l1b': ObjectCapture(rf.Level1b),
            'retrieval_components': ObjectCapture(creator.retrieval.RetrievalComponents),
        }
        config_def.update(captured_objects)
        config_def['order'] += list(captured_objects.keys())

        return config_def

    def config_instance(self, **strategy_keywords):

        logger.debug("Loading configration")

        config_def = self.config_definition(**strategy_keywords)
        config_inst = process_config(config_def)

        logger.debug("Configuration processing complete")

        return config_inst

    def attach_logging(self, config_inst):

        iter_log = rf.SolverIterationLog(config_inst['state_vector'])

        if config_inst['solver'] is not None:
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

    def run_forward_model(self, config_inst, step_index=None):

        file_inst = config_inst['file_config']
        if config_inst['state_vector'] is not None and 'retrieval' in file_inst and 'initial_guess' in file_inst['retrieval']:
            # Necessary for jacobians to work. Normally this step is done in the solver creator
            # Do this before attaching output so that state vector can be output correctly
            logger.debug("Enabling jacobians for forward model simulation")
            config_inst['state_vector'].update_state(file_inst['retrieval']['initial_guess'])

        self.attach_output(config_inst, step_index)

        config_inst['forward_model'].radiance_all()


    def update_covariance(self, config_inst):

        logger.debug("Updating covariances for next step using a posteriori covariance")

        if config_inst['retrieval_components'] is None:
            logger.error("Cannot update covariances: retrieval_components not captured")
            return

        retrieval_components = config_inst['retrieval_components']
    
        if config_inst['state_vector'] is None or config_inst['solver'] is None:
            logger.error("Cannot update covariances: solver or state_vector not defined")
            return

        solver = config_inst["solver"]
        state_vector = config_inst["state_vector"]

        if not hasattr(solver, "problem") and not hasattr(solver.problem, "max_a_posteriori"):
            logger.error("Cannot update covariances: Solver is not using MaxAPosteriori problem")
            return

        max_apost = solver.problem.max_a_posteriori

        # Update state vector to push covariances out to components
        state_vector.update_state(state_vector.state, max_apost.a_posteriori_covariance)

        # First make sure we have all the values we need to perform covariance updating
        if len(self.covariance_storage) == 0:
            logger.error("Cannot update covariances: No covariances were saved during the retrieval step. Please added use of the covariance_storage keyword argument to configuration")
            return

        # Pull covariances out of retrieval components back into covariance storage, where they can be
        # used as apriori inputs for the next step
        for rc_name, rc_obj in retrieval_components.items():
            # This is an error indicating there is some disconnect between the names used for covariances and retrieval components
            # They should be the same since covariances are supplied by retrieval component name
            if rc_name not in self.covariance_storage:
                raise Exception("Cannot find covariance for retrieval component {}".format(rc_name))
    
            # Use retrieval flag to subset covariance, this should parallel what is done in CovarianceByComponent to 
            # prepare covariances for input
            if self.covariance_storage[rc_name].shape[0] != rc_obj.statevector_covariance.shape[0]:
                used_indexes = np.ix_(rc_obj.used_flag_value, rc_obj.used_flag_value)
                self.covariance_storage[rc_name][used_indexes] = rc_obj.statevector_covariance[:]
            else:
                self.covariance_storage[rc_name][:] = rc_obj.statevector_covariance[:]

    def retrieval_cleanup(self, config_inst):
        
        # Remove current observers from the state vector since on the next step
        # a difference configuration of the state vector is possible
        config_inst['state_vector'].clear_observers()

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

            if len(step_indexes) > 1:
                # Only require this keyword for multi-step retrievals
                step_keywords['covariance_storage'] = self.covariance_storage

            config_inst = self.config_instance(**step_keywords)

            self.run_solver(config_inst, step_index)

            if len(step_indexes) > 1:
                # No need to update covariance for single step retrievals
                self.update_covariance(config_inst)

            # Clean up from this step before continuing to the next step
            self.retrieval_cleanup(config_inst)
            
            # Force removal from scope to free the allocated memory
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

import os
import sys
import logging
import pprint

import numpy as np
from netCDF4 import Dataset

from refractor_swig import StateMappingAtIndexes

from .config import load_config_module, find_strategy_function
from .configuration_interface import ConfigurationInterface

logger = logging.getLogger(__name__)

class StrategyExecutor(object):
    
    def __init__(self, config_filename=None, config_inst=None,
                 output_filename=None,
                 strategy_filename=None, strategy_list=None,
                 config_args=None):

        if(config_filename):
            logger.debug("Loading configuration from {}".format(config_filename))
            self.config_module = load_config_module(config_filename)
            self._config_inst = None

        else:
            self._config_inst = config_inst
        self.config_args = config_args

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

    def config_instance(self, **strategy_keywords):
        if(self._config_inst):
            logger.debug("Using configuration instance passed to executor")
            return self._config_inst

        logger.debug("Loading configuration")
        config_inst = ConfigurationInterface.create_configuration_instance(self.config_module, config_args=self.config_args, **strategy_keywords)
        logger.debug("Configuration processing complete")

        return config_inst

    def attach_logging(self, config_inst):
        "Instructs the configuration instance to attach logging classes. Present here to be overriden by inheriting classes"

        config_inst.attach_logging()

    def attach_output(self, config_inst, step_index=None, simulation=False):
        "Instructs the configuration instance to attach output classes. Present here to be overriden by inheriting classes"

        config_inst.attach_output(self.output, step_index, simulation)

    def run_solver(self, config_inst, step_index=None):
        self.attach_logging(config_inst)
        self.attach_output(config_inst, step_index)

        if config_inst.solver is None:
            raise Exception("Solver object was not defined in configuration")

        if config_inst.state_vector is None:
            raise Exception("StateVector object was not defined in configuration")
 
        config_inst.solver.solve()

    def run_forward_model(self, config_inst, step_index=None):
        config_inst.set_initial_guess()

        self.attach_logging(config_inst)
        self.attach_output(config_inst, step_index, simulation=True)

        config_inst.radiance_all()

    def retrieval_indexes(self, rc_obj):

        # Replicate old behavior of having retrieval flags for certain state
        # vector parameters where we want to subselect the covariance in case
        # we turn some of the parameters dynamically during a retrieval but
        # supply the whole covariance to the config
        if isinstance(rc_obj.state_mapping, StateMappingAtIndexes):
            return rc_obj.state_mapping.retrieval_indexes
        else:
            return None

    def update_covariance(self, config_inst):
        logger.debug("Updating covariances for next step using a posteriori covariance")
        if config_inst.retrieval_components is None:
            logger.error("Cannot update covariances: retrieval_components not captured")
            return

        retrieval_components = config_inst.retrieval_components
    
        if config_inst.state_vector is None or config_inst.solver is None:
            logger.error("Cannot update covariances: solver or state_vector not defined")
            return

        solver = config_inst.solver
        state_vector = config_inst.state_vector

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
            retrieval_indexes = self.retrieval_indexes(rc_obj)

            if retrieval_indexes is not None:
                used_indexes = np.ix_(retrieval_indexes, retrieval_indexes)
                self.covariance_storage[rc_name][used_indexes] = rc_obj.statevector_covariance[:, :]
            else:
                self.covariance_storage[rc_name][:] = rc_obj.statevector_covariance[:]

    def retrieval_cleanup(self, config_inst):
        # Remove current observers from the state vector since on the next step
        # a difference configuration of the state vector is possible
        config_inst.state_vector.clear_observers()

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

            if config_inst.state_vector is not None:
                logger.debug("Initial State Vector:")
                logger.debug("---------------------")
                logger.debug(f"{config_inst.state_vector}")

            self.run_solver(config_inst, step_index)

            if len(step_indexes) > 1:
                # No need to update covariance for single step retrievals
                self.update_covariance(config_inst)

            config_inst.retrieval_step_completed(step_index)
            if(step_index == len(self.strategy_list) - 1):
                config_inst.retrieval_completed()

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

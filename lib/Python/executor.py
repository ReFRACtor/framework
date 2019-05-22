import os
import sys
import logging
import pprint

from refractor.config import load_config_module, find_config_function, find_strategy_function
from refractor.factory import process_config, creator
from refractor import framework as rf

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
    
    def __init__(self, config_filename, output_filename=None, strategy_filename=None):

        logger.debug("Loading configuration from {}".format(config_filename))
        self.config_module = load_config_module(config_filename)

        if strategy_filename:
            logger.debug("Loading strategy from {}".format(strategy_filename))
            self.strategy_module = load_config_module(strategy_filename)
        else:
            self.strategy_module = None

    def config_definition(self, **strategy_keywords):

        config_func = find_config_function(self.config_module)

        # Augment loaded configuration to help capture objects required for ouput
        config_def = {
            'creator': creator.base.ParamPassThru, 
            'order': ['file_config', 'forward_model', 'solver', 'state_vector'],
            'file_config': config_func(**strategy_keywords),

            # Capture certain objects without depending on the structure of the configuration
            'forward_model': { 'creator': ObjectCapture(rf.ForwardModel) },   
            'solver': { 'creator': ObjectCapture(rf.IterativeSolver) },   
            'state_vector': { 'creator': ObjectCapture(rf.StateVector) },   
        }

        return config_def

    def strategy_list(self):

        if self.strategy_module is not None:
            return find_strategy_function(self.strategy_module)()
        else:
            # A single empty strategy
            return [ {} ]

    def run_solver(self, config_inst):
        
        config_inst['solver'].solve()
        config_inst['state_vector'].clear_observers()

    def run_forward_model(self, config_inst):

        config_inst['forward_model'].radiance_all()

    def execute_retrieval(self, step_indexes=None):

        strategy = self.strategy_list()
        logger.debug("{} strategy steps loaded".format(len(strategy)))

        if step_indexes is None or len(step_indexes) == 0:
            step_indexes = range(len(strategy))
        else:
            logger.debug("Using {} steps".format(len(step_indexes)))

        for step_index in step_indexes:
            if step_index < 0 or step_index >= len(strategy):
                raise Exception("Invalid step index: {}".format(step_index))

            logger.info("Executing retrieval step {}".format(step_index))

            step_keywords = strategy[step_index]

            logger.debug("Strategy step options:")
            logger.debug(pprint.pformat(step_keywords, indent=2))

            config_def = self.config_definition(**step_keywords)
            config_inst = process_config(config_def)

            self.run_solver(config_inst)

            del config_inst

    def execute_simulation(self, step_index=0):

        strategy = self.strategy_list()

        if step_index < 0 or step_index >= len(strategy):
            raise Exception("Invalid step index: {}".format(step_index))

        logger.info("Executing simulation for step {}".format(step_index))

        step_keywords = strategy[step_index]
        
        config_def = self.config_definition(**step_keywords)
        config_inst = process_config(config_def)

        self.run_forward_model(config_inst)

        del config_inst

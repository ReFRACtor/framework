from abc import ABC, abstractmethod
from refractor_swig import (SolverIterationLog, RtAtmosphere, ForwardModel, IterativeSolver,
                            StateVector, Level1b)
from .config import find_config_function
from ..factory import process_config, creator
from ..output.atmosphere import AtmosphereOutput
from ..output.radiance import (ForwardModelRadianceOutput,
                              ObservationRadianceOutput)
from ..output.solver import SolverIterationOutput
from ..output.state_vector import (StateVectorOutputRetrieval,
                                  StateVectorOutputSimulation)


class ConfigurationInterface(ABC):
    '''This class provides the interface needed by a configuration class 
    by StrategyExecutor.  Note that this just specifies the interface 
    needed, we use the standard python duck typing and a class doesn't
    need to actually derive from this one if it isn't convenient. 

    Note that the creator for this should take a list (possibly empty) of
    strategy_keywords'''
    @property
    @abstractmethod
    def state_vector(self):
        '''State vector used by solver.'''
        pass

    @property
    @abstractmethod
    def solver(self):
        '''Solver to use'''
        pass
    
    @abstractmethod
    def set_initial_guess(self):
        '''Set the state vector to the initial guess.'''
        pass

    @abstractmethod
    def radiance_all(self):
        '''Calculate all the radiances. This doesn't really do anything on
        its own, but if we have been attached to output this would produce
        output. This is used by run_forward_model of StrategyExecutor'''
        pass

    def attach_output(self, output, step_index):
        '''Attach output for the given step index. 

        Note this typically uses the framework Output class, which allows for
        decentralized output. This was an advantage for the OCO-2 development,
        where the output depended on the configuration. But for a particular
        instrument/configuration this might be overly complicated. 

        Some configurations may choose to not do any attachment, and instead
        just respond to retrieval_step_completed and/or retrieval_completed'''
        pass

    def retrieval_step_completed(self, step_index):
        '''Called when a strategy step has been completed, so any desired
        output can be written out, or whatever other response to a a completed
        step is desired.'''
        pass

    def retrieval_completed(self):
        '''Called when a retrieval has been completed, so any desired
        output can be written out, or whatever other response to a a completed
        retrieval is desired.'''
        pass

    def attach_logging(self):
        '''Attach logging. The default implementation uses 
        SolverIterationLog to write out information about the state vector
        as we iterate, but a derived class can specify whatever is desired'''
        iter_log = SolverIterationLog(self.state_vector)

        if self.solver is not None:
            self.solver.add_observer_and_keep_reference(iter_log)

    @classmethod
    def create_configuration_instance(cls, config_module, **strategy_keywords):
        '''Return a configuration instance using the given module (generally
        loaded by config.load_config_module function). This either 
        uses a supplied class, or if a function is defined this is used 
        with ConfigurationCreator'''
        # Despite the name, this is either a class or function
        t = find_config_function(config_module)
        if isinstance(t, type):
            return t(**strategy_keywords)
        return ConfigurationCreator(t, **strategy_keywords)

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

    
class ConfigurationCreator(ConfigurationInterface):
    '''Implementation of ConfigurationInterface that uses a dictionary of
    creator function to create all the objects'''
    def __init__(self, config_func=None, config_dict=None, **strategy_keywords):
        '''Create a ConfigurationCreator. You can either pass in a 
        config_func which is called with the given **strategy_keywords,
        or you can pass in a config_dict which is a dictionary of creators 
        (so like the results of calling this function). Depending on the 
        context, either can be useful. '''
        # Augment loaded configuration to help capture objects
        # required for output
        config_def = {
            'creator': creator.base.ParamPassThru, 
            'order': ['file_config'],
            'file_config': config_func(**strategy_keywords) if config_func else config_dict,
        }

        # Capture certain objects without depending on the structure
        # of the configuration
        captured_objects = {
            'atmosphere': ObjectCapture(RtAtmosphere),
            'forward_model': ObjectCapture(ForwardModel),
            'solver': ObjectCapture(IterativeSolver),
            'state_vector': ObjectCapture(StateVector),
            'l1b': ObjectCapture(Level1b),
            'retrieval_components': ObjectCapture(creator.retrieval.RetrievalComponents),
        }
        config_def.update(captured_objects)
        config_def['order'] += list(captured_objects.keys())
        self.config_inst = process_config(config_def)

    @property
    def file_config(self):
        return self.config_inst.file_config

    @property
    def atmosphere(self):
        return self.config_inst.atmosphere
    
    @property
    def forward_model(self):
        return self.config_inst.forward_model

    @property
    def l1b(self):
        return self.config_inst.l1b

    @property
    def retrieval_components(self):
        return self.config_inst.retrieval_components
    
    @property
    def state_vector(self):
        '''State vector used by solver.'''
        return self.config_inst.state_vector

    @property
    def solver(self):
        '''Solver to use'''
        return self.config_inst.solver

    @property
    def initial_guess(self):
        '''Return the initial guess given to the solver'''
        return self.solver.problem.parameters

    @property
    def observation(self):
        '''Return the observation given to the problem'''
        return self.solver.problem.max_a_posteriori.observation[0]
    
    @property
    def a_priori(self):
        '''Return the a_priori given the to the solver problem.'''
        return self.solver.problem.max_a_posteriori.a_priori_params

    @property
    def a_priori_covariance(self):
        '''Return the a_priori covariance given the to the solver problem.'''
        return self.solver.problem.max_a_posteriori.a_priori_cov
        
    def set_initial_guess(self):
        '''Set the state vector to the initial guess.'''
        if self.state_vector is not None:
            self.state_vector.update_state(self.file_config.retrieval.initial_guess)

    def radiance_all(self):
        '''Calculate all the radiances. This doesn't really do anything on
        its own, but if we have been attached to output this would produce
        output. This is used by run_forward_model of StrategyExecutor'''
        return self.forward_model.radiance_all()

    def attach_output(self, output, step_index):
        '''Attach output for the given step index'''
        if output is None:
            return

        rad_out = ForwardModelRadianceOutput(output, step_index, self.solver)
        self.forward_model.add_observer_and_keep_reference(rad_out)

        if(self.l1b is not None and self.solver is not None):
            obs_out = ObservationRadianceOutput(output, step_index, self.l1b, self.forward_model)
            self.solver.add_observer_and_keep_reference(obs_out)

        if self.atmosphere is not None and self.solver is not None:
            atm_out = AtmosphereOutput(output, step_index, self.atmosphere)
            self.solver.add_observer_and_keep_reference(atm_out)

        if self.state_vector is not None and self.solver is not None:
            solver_out = SolverIterationOutput(output, self.solver, step_index)
            self.solver.add_observer_and_keep_reference(solver_out)

            sv_out = StateVectorOutputRetrieval(output, step_index, self.state_vector)
            self.solver.add_observer_and_keep_reference(sv_out)

        elif self.state_vector is not None and self.solver is None:
            sv_out = StateVectorOutputSimulation(output, step_index, self.state_vector)

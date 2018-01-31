from collections import OrderedDict

import numpy as np

from .base import Creator
from .. import param

from refractor import framework as rf

class NLLSRetrieval(Creator):

    retrieval_components = param.Dict()
    state_vector = param.InstanceOf(rf.StateVector)
    initial_guess = param.Array(dims=1)
    a_priori = param.Array(dims=1)
    covariance = param.Array(dims=1)
    solver = param.InstanceOf(rf.NLLSSolver)

    def __init__(self, *vargs, **kwargs):
        super().__init__(*vargs, **kwargs)

    def create(self, **kwargs):

        # The order these are accessed matters, store values into common for subsequent steps
        retrieval_components = self.common_store["retrieval_components"] = self.retrieval_components()
        state_vector = self.common_store["state_vector"] = self.state_vector()
        initial_guess = self.common_store["initial_guess"] = self.initial_guess()
        a_priori = self.common_store["a_priori"] = self.a_priori()


        return {
            'retrieval_components': retrieval_components,
            'state_vector': state_vector,
            'intial_guess': initial_guess,
            'a_priori': a_priori,
            'covariance': None,
            'solver': None,
        }

class SVObserverComponents(Creator):

    def __init__(self, *vargs, **kwargs):
        super().__init__(*vargs, **kwargs)

        self.register_to_receive(rf.SubStateVectorObserver)

        self.retrieval_components = OrderedDict()

    def receive(self, rec_obj):

        if hasattr(rec_obj, "sub_state_identifier") and rec_obj.coefficient.value.shape[0] > 0:
            ss_iden = rec_obj.sub_state_identifier
            self.retrieval_components[ss_iden] = rec_obj
 
    def create(self, **kwargs):
        return self.retrieval_components

class StateVector(Creator):

    retrieval_components = param.Dict()

    def __init__(self, *vargs, **kwargs):
        super().__init__(*vargs, **kwargs)

        self.register_to_receive(rf.StateVectorObserver)

        self.sv_observers = []

    def receive(self, rec_obj):

        if not hasattr(rec_obj, "sub_state_identifier") and rec_obj not in self.sv_observers:
            self.sv_observers.append(rec_obj)

    def create(self, **kwargs):
        sv = rf.StateVector()

        # Register retrieval components first as state vector observers
        for observer in self.retrieval_components().values():
            sv.add_observer(observer)
            
        # Add remaining non retrieval components as observers
        for observer in self.sv_observers:
            sv.add_observer(observer)

        return sv

class InitialGuessFromSV(Creator):

    retrieval_components = param.Dict()
    state_vector = param.InstanceOf(rf.StateVector)

    def create(self, **kwargs):

        sv = self.state_vector()

        # Gather initial guess from retrieval components
        ig_values = []
        for ret_component in self.retrieval_components().values():
            used_indexes = np.nonzero(ret_component.used_flag_value)
            ig_values.append(ret_component.coefficient.value[used_indexes])

        ig = np.concatenate(ig_values)

        if ig.shape[0] != sv.observer_claimed_size:
            raise ValueError("The initial guess vector size %d does not match expected state vector size %d" % (ig.shape[0], sv.observer_claimed_size))
        
        sv.update_state(ig)

        return ig

class AprioriFromIG(Creator):

    initial_guess = param.Array(dims=1)

    def create(self, **kwargs):

        return self.initial_guess()

class CovarianceByComponent(Creator):

    retrieval_components = param.Dict()
    values = param.Dict() 

    def create(self, **kwargs):

        cov_inputs = self.values()

        # Gather covariance arrays and do input checking
        covariances = []
        total_len = 0
        for rc_name in self.retrieval_components.keys():
            if rc_name not in cov_inputs:
                raise param.ParamError("CovarianceByComponent: values argument is missing data for this retrieval component: %s" % rc_name)
            
            rc_cov = cov_inputs[rc_name]

            if not hasattr("shape", rc_cov):
                raise param.ParamError("CovarianceByComponent: value for retrieval component must be a numpy array: %s" % rc_name)

            if len(rc_cov.shape) != 2:
                raise param.ParamError("CovarianceByComponent: shape of array for retrieval component must be 2 dimensional: %s" % rc_name)

            if rc_cov.shape[0] != rc_cov.shape[1]:
                raise param.ParamError("CovarianceByComponent: array for retrieval component must be a square matrix: %s" % rc_name)

            total_len += rc_cov.shape[0]
            covariances.append(rc_cov)

        # Concatenate covariances along the diagonal
        total_covariance = np.zeros((total_len, total_len), dtype=float)

        offset = 0
        for idx, cov in enumerate(covariances):
            total_covariance[offset:cov.shape[0], offset:cov.shape[1]] = cov
            offset += cov.shape[0]

        return total_covariance

class NLLSMaxAPosteriori(Creator):

    l1b = param.InstanceOf(rf.Level1b)
    instrument = param.InstanceOf(rf.Instrument)
    forward_model = param.InstanceOf(rf.ForwardModel)

    state_vector = param.InstanceOf(rf.StateVector)

    a_priori = param.Array(dims=1)
    covariance = param.Array(dims=2)

    max_cost_function_calls = param.Scalar(int)
    dx_tol_abs = param.Scalar(float)
    dx_tol_rel = param.Scalar(float)
    g_tol_abs = param.Scalar(float)

    def create(self, **kwargs):

        fm = self.forward_model()
        observation = rf.ObservationLevel1b(self.l1b(), self.instrument(), fm.spectral_grid)

        stat_method = rf.MaxAPosterioriStandard(fm, observation, self.state_vector(), self.a_priori(), self.covariance())

        opt_problem = rf.NLLSMaxAPosteriori(stat_method, True)

        solver = NLLSSolverGSLLMSDER(self.max_cost_function_calls(), 
                                     self.dx_tol_abs(), self.dx_tol_rel(), self.g_tol_abs(),
                                     opt_problem, True)

        return solver

import re
from collections import OrderedDict

import numpy as np

from .base import Creator
from .. import param

from refractor import framework as rf

# Class allow capturing the list of retrieval components from the configuration system elsewhere
class RetrievalComponents(OrderedDict):
    pass

class RetrievalBaseCreator(Creator):

    retrieval_components = param.Dict()
    state_vector = param.InstanceOf(rf.StateVector)
    initial_guess = param.Array(dims=1)

    def __init__(self, *vargs, **kwargs):
        super().__init__(*vargs, **kwargs)

    def create(self, **kwargs):

        # The order these are accessed matters, store values into common for subsequent steps
        retrieval_components = self.common_store["retrieval_components"] = self.retrieval_components()
        state_vector = self.common_store["state_vector"] = self.state_vector()
        initial_guess = self.common_store["initial_guess"] = self.initial_guess()

        return {
            'retrieval_components': retrieval_components,
            'state_vector': state_vector,
            'initial_guess': initial_guess,
        }
  
class NLLSRetrieval(RetrievalBaseCreator):

    a_priori = param.Array(dims=1)
    covariance = param.Array(dims=2)
    solver = param.Choice(param.InstanceOf(rf.NLLSSolver), param.InstanceOf(rf.ConnorSolver))

    def __init__(self, *vargs, **kwargs):
        super().__init__(*vargs, **kwargs)

    def create(self, **kwargs):
        retrieval_objects = super().create(**kwargs)

        # The order these are accessed matters, store values into common for subsequent steps
        a_priori = self.common_store["a_priori"] = self.a_priori()
        covariance = self.common_store["covariance"] = self.covariance()
        solver = self.common_store["solver"] = self.solver()

        retrieval_objects.update({
            'a_priori': a_priori,
            'covariance': covariance,
            'solver': solver,
        })

        return retrieval_objects

class SVObserverComponents(Creator):

    exclude = param.Iterable(default=None, required=False)
    include = param.Iterable(default=None, required=False)
    order = param.Iterable(default=[], required=False)

    # Alternative names for state vector components:
    # alt_name -> sub_state_identifier
    alt_component_names = {
        "TSUR": "surface_temperature",
        "TATM": "^temperature",
        "EMIS": "emissivity",
    }

    def __init__(self, *vargs, **kwargs):
        super().__init__(*vargs, **kwargs)

        self.register_to_receive(rf.SubStateVectorObserver)

        self.retrieval_components = OrderedDict()

    def receive(self, rec_obj):

        if hasattr(rec_obj, "sub_state_identifier") and rec_obj.coefficient.value.shape[0] > 0:
            ss_iden = rec_obj.sub_state_identifier
            self.retrieval_components[ss_iden] = rec_obj

    def in_component_list(self, iden, comp_list):

        for comp_re in comp_list:
            if comp_re in self.alt_component_names:
                # Convert alternative names to one that will match against the sub_state_identifier strings
                comp_re = self.alt_component_names[comp_re]
            elif not re.search('/', comp_re) and comp_re in self.sub_comp_names:
                # If the items does not contain a slash but matches something after a slash in the list of 
                # possible components, prepend a slash to prevent spurious matches, ie, "O3" matching "HNO3"
                comp_re = "/" + comp_re

            if re.search(comp_re, iden):
                return True
        return False

    def is_included(self, iden):

        include_list = self.include()

        if include_list is None:
            return True
        else:
            return self.in_component_list(iden, include_list)

    def is_excluded(self, iden):

        exclude_list = self.exclude()

        if exclude_list is None:
            return False
        else:
            return self.in_component_list(iden, exclude_list)

    def create(self, **kwargs):
        
        # Names of retrieval components after a slash, for example, "CO2" after "absorber_scaled/CO2"
        # Used to disambiguate matches
        self.sub_comp_names = [ c.split("/")[-1] for c in filter(lambda c: re.search('/', c), self.retrieval_components.keys()) ]
        
        filtered_components = []
        for iden, obj in self.retrieval_components.items():
            if self.is_included(iden) and not self.is_excluded(iden):
                filtered_components.append( (iden, obj) )

        # Sort by the items listed in the order parameter
        # falling back to the order of items originally recieved
        order = self.order()
        def sort_key(comp_tuple):
            found_order = None
            for idx, ord_item in enumerate(order):
                if re.search(ord_item, comp_tuple[0]):
                    found_order = idx
                    break

            if found_order is not None:
                return found_order
            else:
                # Add length of order array so items without a specific order
                # come after any that match the specific ordering
                dflt_order = len(order) + filtered_components.index(comp_tuple)
                return dflt_order

        return RetrievalComponents(sorted(filtered_components, key=sort_key))

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

        if len(ig_values) == 0:
            raise param.ParamError("InitialGuessFromSV: No initial guess values available as identified by the retrieval components")

        ig = np.concatenate(ig_values)

        if ig.shape[0] != sv.observer_claimed_size:
            raise ValueError("The initial guess vector size %d does not match expected state vector size %d" % (ig.shape[0], sv.observer_claimed_size))

        return ig

class AprioriFromIG(Creator):

    initial_guess = param.Array(dims=1)

    def create(self, **kwargs):

        return self.initial_guess()

class CovarianceByComponent(Creator):

    retrieval_components = param.Dict()
    values = param.Dict() 
    interstep_storage = param.Dict(required=False)

    def create(self, **kwargs):

        storage = self.interstep_storage()
        if storage is not None:
            if len(storage) == 0:
                storage.update(self.values())
            cov_inputs = storage
        else:
            cov_inputs = self.values()

        # Gather covariance matrices and do input checking
        covariances = []
        total_len = 0
        for rc_name, rc_obj in self.retrieval_components().items():
            if rc_name not in cov_inputs:
                raise param.ParamError("CovarianceByComponent: covariance values argument is missing data for this retrieval component: %s" % rc_name)
            
            rc_cov = cov_inputs[rc_name]

            if not hasattr(rc_cov, "shape"):
                raise param.ParamError("CovarianceByComponent: value for retrieval component covariance must be a numpy array: %s" % rc_name)

            if len(rc_cov.shape) != 2:
                raise param.ParamError("CovarianceByComponent: shape of array for retrieval component covariance must be 2: %s" % rc_name)

            if rc_cov.shape[0] != rc_cov.shape[1]:
                raise param.ParamError("CovarianceByComponent: array for retrieval component must be a square matrix: %s" % rc_name)

            # Subset a larger covariances if it is the same size as the flags
            flag = rc_obj.used_flag_value
            if flag.shape[0] == rc_cov.shape[0]:
                used_indexes = np.nonzero(flag)
                used_cov = rc_cov[np.ix_(used_indexes[0], used_indexes[0])]
            else:
                used_cov = rc_cov

            # Check that the matrix is postive definite to help reduce debugging time by highlighting this fact before
            # the full Cholesky matrix decomposition within the framework fails
            if not np.all(np.linalg.eigvals(used_cov) > 0):
                raise param.ParamError("CovarianceByComponent: covariance for {} is not positive definite".format(rc_name))

            total_len += used_cov.shape[0]
            covariances.append(used_cov)

        # Concatenate covariances along the diagonal
        total_covariance = np.zeros((total_len, total_len), dtype=float)

        offset = 0
        for idx, cov in enumerate(covariances):
            total_covariance[offset:offset+cov.shape[0], offset:offset+cov.shape[1]] = cov
            offset += cov.shape[0]

        return total_covariance

class LegacyConnorSolver(Creator):

    l1b = param.InstanceOf(rf.Level1b)
    instrument = param.InstanceOf(rf.Instrument)
    forward_model = param.InstanceOf(rf.ForwardModel)

    state_vector = param.InstanceOf(rf.StateVector)

    initial_guess = param.Array(dims=1)
    a_priori = param.Array(dims=1)
    covariance = param.Array(dims=2)

    max_iteration = param.Scalar(int)
    max_divergence = param.Scalar(int)
    max_chisq = param.Scalar(float)
    threshold = param.Scalar(float)
    gamma_initial = param.Scalar(float)

    def create(self, **kwargs):
        sv = self.state_vector()
        fm = self.forward_model()

        observation = rf.ObservationLevel1b(self.l1b(), self.instrument(), fm.spectral_grid)
    
        cost_func = rf.ConnorCostFunction(sv, fm, observation)
        conv = rf.ConnorConvergence(fm,
                                 self.threshold(), 
                                 self.max_iteration(), 
                                 self.max_divergence(), 
                                 self.max_chisq())

        # Create wrapper function to make the ConnorSolver interface match that of IterativeSolver
        # for at least the solve routine
        def create_solve_wrapper(conn_solver, conn_creator):
            orig_solve_func = conn_solver.solve

            def solve_wrapper(**vargs):
                "Automatically pass initial guess, apriori and covariance to solve routine if no arguments supplied"

                if len(vargs) == 3:
                    return orig_solve_func()
                elif len(vargs) == 0:
                    return orig_solve_func(conn_creator.initial_guess(), conn_creator.a_priori(), conn_creator.covariance())
                else:
                    raise Exception("Wrong number of arguments to ConnorSolver::solve. Expected 0 or 3 arguments.")

            return solve_wrapper

        solver = rf.ConnorSolver(cost_func, conv, self.gamma_initial())
        solver.solve = create_solve_wrapper(solver, self)

        # Use add_observer_and_keep_reference to keep observer from going away as soon
        # as this function ends, since observers are normally kept by a weak reference
        iter_log = rf.ConnorIterationLog(sv)
        solver.add_observer_and_keep_reference(iter_log)

        return solver

class MaxAPosterioriBase(Creator):

    l1b = param.InstanceOf(rf.Level1b)
    instrument = param.InstanceOf(rf.Instrument)
    forward_model = param.InstanceOf(rf.ForwardModel)

    state_vector = param.InstanceOf(rf.StateVector)

    initial_guess = param.Array(dims=1)
    a_priori = param.Array(dims=1)
    covariance = param.Array(dims=2)

    max_cost_function_calls = param.Scalar(int)

    def init_state_vector(self):

        self.state_vector().update_state(self.a_priori(), self.covariance())

    def opt_problem(self):
        fm = self.forward_model()

        observation = rf.ObservationLevel1b(self.l1b(), self.instrument(), fm.spectral_grid)

        apriori = self.a_priori()
        cov = self.covariance()

        if apriori.shape[0] != cov.shape[0] or apriori.shape[0] != cov.shape[1]:
            raise param.ParamError(f"The number of rows and columns of a priori covariance matrix: {cov.shape} do not equal the size of the a priori parameter values array: {apriori.shape[0]}")

        stat_method = rf.MaxAPosterioriStandard(fm, observation, self.state_vector(), apriori, cov)

        opt_problem = rf.NLLSMaxAPosteriori(stat_method)

        # Initialize solver intitial guess
        opt_problem.parameters = self.initial_guess()

        return opt_problem

class NLLSSolverGSLLMSDER(MaxAPosterioriBase):
    """GSLLM<with S>DER solver

    J. J. More's version of the Levenberg-Marquardt NLLS solver with a diagonal weighting matrix"""
 

    dx_tol_abs = param.Scalar(float)
    dx_tol_rel = param.Scalar(float)
    g_tol_abs = param.Scalar(float)

    def create(self, **kwargs):

        solver = rf.NLLSSolverGSLLMSDER(self.opt_problem(),
                self.max_cost_function_calls(), self.dx_tol_abs(), self.dx_tol_rel(), self.g_tol_abs(),
                False)

        self.init_state_vector()

        return solver

class NLLSSolverGSLLMDER(MaxAPosterioriBase):
    """GSLLM<no S>DER solver.

    J. J. More's version of the Levenberg-Marquardt NLLS solver with one difference.  
    The diagonal weighing matrix is replaced with the identity matrix.
    """

    dx_tol_abs = param.Scalar(float)
    dx_tol_rel = param.Scalar(float)
    g_tol_abs = param.Scalar(float)

    def create(self, **kwargs):

        solver = rf.NLLSSolverGSLLMDER(self.opt_problem(),
                self.max_cost_function_calls(), self.dx_tol_abs(), self.dx_tol_rel(), self.g_tol_abs(),
                False)

        self.init_state_vector()

        return solver

class ConnorSolverMAP(MaxAPosterioriBase):

    threshold = param.Scalar(float)
    max_iteration = param.Scalar(int)
    max_divergence = param.Scalar(int)
    max_chisq = param.Scalar(float)
    gamma_initial = param.Scalar(float)

    def create(self, **kwargs):

        conv = rf.ConnorConvergence(self.forward_model(), self.threshold(), self.max_iteration(), self.max_divergence(), self.max_chisq())

        solver = rf.ConnorSolverMAP(self.opt_problem(), conv,
                self.max_cost_function_calls(), True,
                self.gamma_initial())

        self.init_state_vector()

        return solver
 
class NLLSSolverLM(MaxAPosterioriBase):

    max_iteration = param.Scalar(int)

    dx_tol_abs = param.Scalar(float, default=1e-6)
    dx_tol_rel = param.Scalar(float, default=1e-6)
    g_tol_abs = param.Scalar(float, default=1e-6)
    g_tol_rel = param.Scalar(float, default=1e-6)

    min_W = param.Scalar(float, required=False)
    tr_rad_tol = param.Scalar(float, required=False)
    tr_rad = param.Scalar(float, required=False)
    cr_ratio_tol = param.Scalar(float, required=False)

    def create(self, **kwargs):

        # Class comes with defaults, only overwrite the defaults if the config has a value supplied
        opts = rf.NLLSSolverLMOptions()

        if self.min_W() is not None:
            opts.min_W = self.min_W()

        if self.tr_rad_tol() is not None:
            opts.tr_rad_tol = self.tr_rad_tol()

        if self.tr_rad() is not None:
            opts.tr_rad = self.tr_rad()

        if self.cr_ratio_tol() is not None:
            opts.cr_ratio_tol = self.cr_ratio_tol()

        solver = rf.NLLSSolverLM(self.opt_problem(), self.max_iteration(), opts, self.dx_tol_abs(), self.dx_tol_rel(), self.g_tol_abs(), self.g_tol_rel())

        self.init_state_vector()

        return solver

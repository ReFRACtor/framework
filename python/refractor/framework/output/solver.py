import logging

import numpy as np
from netCDF4 import stringtochar

import refractor.framework as rf
from .base import OutputBase

logger = logging.getLogger(__name__)

class SolverIterationOutput(rf.ObserverIterativeSolver, OutputBase):
    "Notifies output classes the current retrieval iteration number when it changes"

    def __init__(self, output, config_solver, step_index):
        # Required to initialize director
        rf.ObserverIterativeSolver.__init__(self)

        self.output = output
        self.step_index = step_index
        self.solver = config_solver

    def notify_update(self, iter_solver):

        # We do not use the iter_solver passed when an update
        # occurs because it will always be a IterSolver object
        # Instead use the solver object passed in at object
        # construction so we can use solver implementation
        # features (after checking for their presence)

        iter_index = self.solver.num_accepted_steps

        # Record number of steps encountered
        base_group_name = self.base_group_name(iter_index)
        base_group = self.output.createGroup(base_group_name)

        if "number_steps" in base_group.variables:
            num_steps = base_group["number_steps"]
        else:
            num_steps = base_group.createVariable("number_steps", int)

        num_steps[...] = self.step_index+1
        
        # Record iteration index at the root of each step to keep a record
        # of the number of iterations encountered

        step_group_name = self.step_group_name(self.step_index, iter_index)
        step_root_group = self.output.createGroup(step_group_name)

        if "number_iterations" in step_root_group.variables:
            num_iters = step_root_group["number_iterations"]
        else:
            num_iters = step_root_group.createVariable("number_iterations", int)

        num_iters[...] = iter_index

        # Solver output is per iteration
        iter_group_name = self.iter_step_group_name(self.step_index, iter_index)

        status_dim = "status_string_s{}".format(self.step_index+1)
        if status_dim not in self.output.dimensions:
            self.output.createDimension(status_dim)

        # Solver output
        solver_group_name = iter_group_name + "Solver"
        solver_group = self.output.createGroup(solver_group_name)

        logger.debug("Storing solver information into output file at {}".format(solver_group_name))
 
        if "cost_function" in solver_group.variables:
            cost_function = solver_group["cost_function"]
        else:
            cost_function = solver_group.createVariable("cost_function", float)
        cost_function[...] = self.solver.cost_at_accepted_points[iter_index]

        if "status" in solver_group.variables:
            status = solver_group["status"]
        else:
            status = solver_group.createVariable("status", "S1", (status_dim))
        status_len = len(self.solver.status_str)
        status[:] = stringtochar(np.array(self.solver.status_str, f"S{status_len}"))

        if "num_accepted" in solver_group.variables:
            num_accepted = solver_group["num_accepted"]
        else:
            num_accepted = solver_group.createVariable("iteration_index", int)
        num_accepted[...] = self.solver.num_accepted_steps

        if hasattr(self.solver, "problem") and hasattr(self.solver.problem, "max_a_posteriori"):
            apost_cov = self.solver.problem.max_a_posteriori.a_posteriori_covariance

            # This should be the same as in StateVectorOutput, since they are sized the same
            sv_dim = "state_vector_s{}".format(self.step_index+1)
            if sv_dim not in self.output.dimensions:
                self.output.createDimension(sv_dim, apost_cov.shape[0])

            if "covariance_a_posteriori" in solver_group.variables:
                cov = solver_group["covariance_a_posteriori"]
            else:
                cov = solver_group.createVariable("covariance_a_posteriori", float, (sv_dim, sv_dim))
            cov[:, :] = apost_cov

        if hasattr(self.solver, "gradient_at_accepted_points") and len(self.solver.gradient_at_accepted_points) > 0:
            grad_at_iter = self.solver.gradient_at_accepted_points[iter_index-1]

            # This should be the same as in StateVectorOutput, since they are sized the same
            sv_dim = "state_vector_s{}".format(self.step_index+1)
            if sv_dim not in self.output.dimensions:
                self.output.createDimension(sv_dim, grad_at_iter.shape[0])

            if "gradient" in solver_group.variables:
                gradient = solver_group["gradient"]
            else:
                gradient = solver_group.createVariable("gradient", float, (sv_dim,))

            gradient[...] = grad_at_iter

import logging

import numpy as np
from netCDF4 import stringtochar

from refractor import framework as rf
from .base import OutputBase

logger = logging.getLogger(__name__)

class SolverIterationOutput(rf.ObserverIterativeSolver, OutputBase):
    "Notifies output classes the current retrieval iteration number when it changes"

    dimension_names = {
        "state_vector": "state_vector",
        "sv_string": "sv_string",
        "status_string": "status_string",
    }

    def __init__(self, output, step_index, state_vector):
        # Required to initialize director
        rf.ObserverIterativeSolver.__init__(self)

        self.output = output
        self.step_index = step_index

        self.state_vector = state_vector

        self.create_dimensions(self.step_index)

    def notify_update(self, solver):

        iter_index = solver.num_accepted_steps
        base_group_name = self.iter_step_group_name(self.step_index, iter_index)

        sv_dim = self.dimensions["state_vector"].name
        str_dim = self.dimensions["sv_string"].name
        status_dim = self.dimensions["status_string"].name

        # StateVector output
        sv_group_name = base_group_name + "StateVector"
        sv_group = self.output.createGroup(sv_group_name)

        logger.debug("Storing state vector information into output file at {}".format(sv_group_name))
 
        if "names" in sv_group.variables:
            sv_names = sv_group["names"]
        else:
            sv_names = sv_group.createVariable("names", "S1", (sv_dim, str_dim))
        sv_names[:] = stringtochar(np.array(self.state_vector.state_vector_name, dtype="S3"))

        if "values" in sv_group.variables:
            sv_values = sv_group["values"]
        else:
            sv_values = sv_group.createVariable("values", float, (sv_dim))
        sv_values[:] = self.state_vector.state

        # Solver output
        solver_group_name = base_group_name + "Solver"
        solver_group = self.output.createGroup(solver_group_name)

        logger.debug("Storing solver information into output file at {}".format(solver_group_name))
 
        if "accepted_points" in solver_group.variables:
            accepted_points = solver_group["accepted_points"]
        else:
            accepted_points = solver_group.createVariable("accepted_points", float, (sv_dim))
        accepted_points[...] = solver.accepted_points[iter_index]

        if "cost_function" in solver_group.variables:
            cost_function = solver_group["cost_function"]
        else:
            cost_function = solver_group.createVariable("cost_function", float)
        cost_function[...] = solver.cost_at_accepted_points[iter_index]

        if "status" in solver_group.variables:
            status = solver_group["status"]
        else:
            status = solver_group.createVariable("status", "S1", (status_dim))
        status[:] = stringtochar(np.array(solver.status_str, "S3"))

        if "num_accepted" in solver_group.variables:
            num_accepted = solver_group["num_accepted"]
        else:
            num_accepted = solver_group.createVariable("iteration_index", int)
        num_accepted[...] = solver.num_accepted_steps

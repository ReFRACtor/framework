import logging

import numpy as np

from netCDF4 import stringtochar

logger = logging.getLogger(__name__)

class OutputBase(object):

    def __init__(self, output):
        self.output = output

    def base_group_name(self, iter_index=None, prefix=""):
        # Determine where to store values
        group_name = prefix

        if len(prefix) > 0 and not prefix[-1] == "/":
            group_name += "/"

        # Without an iteration index assume output is for a simulation
        if iter_index is not None:
            group_name += "Retrieval/"
        else:
            group_name += "Simulation/"

        return group_name

    def step_group_name(self, step_index=None, iter_index=None, prefix=""):

        group_name = self.base_group_name(iter_index, prefix)

        if step_index is not None:
            group_name += "Step_{}/".format(step_index + 1)

        return group_name

    def iter_step_group_name(self, step_index=None, iter_index=None, prefix=""):

        group_name = self.step_group_name(step_index, iter_index, prefix)

        # See if we can determine and iteration index for retrievals
        if iter_index is not None:
            if iter_index == 0:
                group_name += "Initial/"
            else:
                group_name += "Iteration_{}/".format(iter_index)

        return group_name

    def create_dimension(self, dim_name, size=None, step_index=None):

        if step_index is not None:
            dim_name += "_s{}".format(step_index)

        if dim_name not in self.output.dimensions:
            self.output.createDimension(dim_name, size)

        return dim_name

    def create_group(self, grp_name, base_group=None):

        if base_group is None:
            base_group = self.output

        group_obj = base_group.createGroup(grp_name)

        return group_obj

    def create_variable(self, var_name, group, var_type, dims):
    
        if var_name in group.variables:
            var_obj = group[var_name]
        else:
            var_obj = group.createVariable(var_name, var_type, dims)

        return var_obj

    def set_string_variable(self, var_name, group, dims, string_list):

        max_str_len = np.max([ len(n) for n in string_list ])
        var_obj = self.create_variable(var_name, group, "S1", dims)  
        var_obj[:] = stringtochar(np.array(string_list, dtype="S{}".format(max_str_len)))

        return var_obj

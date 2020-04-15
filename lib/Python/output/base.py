import logging

logger = logging.getLogger(__name__)

class OutputBase(object):

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

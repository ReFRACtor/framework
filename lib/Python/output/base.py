import logging

logger = logging.getLogger(__name__)

class OutputBase(object):

    dimension_names = {
        # "internal_name" : "output_name"
    }

    def create_dimensions(self):

        if not hasattr(self, "dimension_names"):
            raise Exception("Output class {} does not have dimension names defined".format(self.__class__.__name__))

        self.dimensions = {}
        for spec_type, dim_name in self.dimension_names.items():
            if dim_name in self.output.dimensions:
                self.dimensions[spec_type] = self.output.dimensions[dim_name]
            else:
                self.dimensions[spec_type] = self.output.createDimension(dim_name, 0)

    def iter_step_group_name(self, step_index=None, iter_index=None, prefix=""):
        # Determine where to store values
        group_name = prefix

        if len(prefix) > 0 and not prefix[-1] == "/":
            group_name += "/"

        # Without an iteration index assume output is for a simulation
        if iter_index is not None:
            group_name += "Retrieval/"
        else:
            group_name += "Simulation/"

        if step_index is not None:
            group_name += "Step_{}/".format(self.step_index)

        # See if we can determine and iteration index for retrievals
        if iter_index is not None:
            if iter_index == 0:
                group_name += "Initial/"
            else:
                group_name += "Iteration_{}/".format(iter_index)

        return group_name

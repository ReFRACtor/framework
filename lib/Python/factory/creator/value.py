from .base import Creator
from .. import param

from refractor import framework as rf

class ArrayWithUnit(Creator):
    "Create an array with unit class from a numpy array and type string"

    value = param.Array()
    units = param.Scalar(str)

    def create(self, **kwargs):

        value = self.value()
        units = self.units()

        num_dims = len(value.shape)
        if num_dims == 3:
            return rf.ArrayWithUnit_double_3(value, units)
        elif num_dims == 2:
            return rf.ArrayWithUnit_double_2(value, units)
        elif num_dims == 1:
            return rf.ArrayWithUnit_double_1(value, units)
        else:
            raise param.ParamError("Unsupported number of dimensions %s for array" % (num_dims))


import numpy as np
import h5py

from .base import Creator
from .. import param

from refractor import framework as rf

class CreatorFlaggedValue(Creator):

    value = param.Choice(param.Array(dims=1), param.ArrayWithUnit(dims=1))
    retrieved = param.Scalar(bool, required=False)

    def retrieval_flag(self):
        val = self.value()
        retrieved = self.retrieved()

        if hasattr(val, "value"):
            # ArrayWithUnit
            val_shape = val.value.shape
        else:
            val_shape = val.shape

        if retrieved is None or retrieved:
            return np.ones(val_shape, dtype=bool)
        else:
            return np.zeros(val_shape, dtype=bool)

class CreatorFlaggedValueMultiChannel(Creator):

    value = param.Choice(param.Array(dims=2), param.ArrayWithUnit(dims=2))
    retrieved = param.Iterable(required=False)

    def retrieval_flag(self):
        val = self.value()
        retrieved = self.retrieved()

        if hasattr(val, "value"):
            # ArrayWithUnit
            val_shape = val.value.shape
        else:
            val_shape = val.shape

        if retrieved is None:
            return np.ones(val_shape, dtype=bool)
        else:
            flags = np.empty(val_shape, dtype=bool)
            for chan_idx, chan_flag in enumerate(retrieved):
                flags[chan_idx, :] = chan_flag
            return flags

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

class LoadValuesFromHDF(Creator):

    filename = param.Scalar(str)

    def create(self, **kwargs):

        values = {}

        contents = h5py.File(self.filename(), "r")

        def visit_datasets(name, obj):
            if type(obj) is h5py.Dataset:
                values[name] = obj[:]

        contents.visititems(visit_datasets)

        return values

import numpy as np
import netCDF4

from .base import Creator
from .. import param

import refractor.framework as rf

class ArrayWithUnit(Creator):
    "Create an array with unit class from a numpy array and type string"

    value = param.Array()
    units = param.Scalar(str)

    def create(self, **kwargs):

        import warnings
        warnings.warn("Instead of using this creator use framework.ArrayWithUnit factory class", DeprecationWarning)

        value = self.value(**kwargs)
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

        contents = netCDF4.Dataset(self.filename(), "r")

        def extract_datasets(group, values):
            for v_name, v_ds in group.variables.items():
                v_path = (group.path + "/" + v_ds.name).lstrip("/")
                values[v_path] = v_ds[:]
            for g_name, child_group in group.groups.items():
                extract_datasets(child_group, values)

        extract_datasets(contents, values)

        return values


class NamedCommonValue(Creator):

    name = param.Scalar(str)

    def create(self, **kwargs):
        self.register_parameter(self.name(), param.AnyValue())
        return self.param(self.name())

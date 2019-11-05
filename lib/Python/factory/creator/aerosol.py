
import os

import numpy as np

from .base import Creator, ParamPassThru
from .value import CreatorFlaggedValue
from .. import param

from refractor import framework as rf


class AerosolDefinition(ParamPassThru):
    "Defines the interface expected for aerosol config defnition blocks, values are pass through as a dictionary"

    extinction = param.InstanceOf(rf.AerosolExtinction)
    properties = param.InstanceOf(rf.AerosolProperty)


class AerosolOptical(Creator):
    "Creates an AbsorberAbsco object that statisfies the AtmosphereCreato;rs absorber value"

    aerosols = param.Iterable()
    pressure = param.InstanceOf(rf.Pressure)
    relative_humidity = param.InstanceOf(rf.RelativeHumidity)
    reference_wn = param.Scalar(int, default=1e4 / 0.755)

    def create(self, **kwargs):

        vec_extinction = rf.vector_aerosol_extinction()
        vec_properties = rf.vector_aerosol_property()

        for aerosol_name in self.aerosols():
            self.register_parameter(aerosol_name, param.Dict())
            aerosol_def = self.param(aerosol_name, aerosol_name=aerosol_name)

            if not "extinction" in aerosol_def:
                raise param.ParamError("exitinction value not in aerosol definition for aerosol: %s" % aerosol_name)

            if not "properties" in aerosol_def:
                raise param.ParamError("peroperties value not in aerosol definition for aerosol: %s" % aerosol_name)

            vec_extinction.push_back(aerosol_def['extinction'])
            vec_properties.push_back(aerosol_def['properties'])

        return rf.AerosolOptical(vec_extinction, vec_properties, self.pressure(), self.relative_humidity(), self.reference_wn())


class AerosolShapeGaussian(CreatorFlaggedValue):

    pressure = param.InstanceOf(rf.Pressure)
    log_space = param.Scalar(bool, default=True)

    def create(self, aerosol_name=None, **kwargs):

        if aerosol_name is None:
            raise param.ParamError("aerosol_name not supplied to creator")

        return rf.AerosolShapeGaussian(self.pressure(), self.retrieval_flag(), self.value(), aerosol_name, not self.log_space())


class AerosolProfileExtinction(CreatorFlaggedValue):

    pressure = param.InstanceOf(rf.Pressure)

    # Callers should specify either log_space or mapping; not both
    log_space = param.Scalar(bool, required=False)
    mapping = param.InstanceOf(rf.Mapping, required=False)

    def create(self, aerosol_name=None, **kwargs):

        if self.log_space() is not None and self.mapping() is not None:
            raise param.ParamError("Specifying both log_space and mapping is ambiguous")
        elif self.log_space() is not None:
            if self.log_space():
                effective_mapping = rf.MappingLog()
            else:
                effective_mapping = rf.MappingLinear()
        elif self.mapping() is not None:
            effective_mapping = self.mapping()
        else:
            effective_mapping = rf.MappingLinear()

        return rf.AerosolExtinctionLevel(self.pressure(), self.retrieval_flag(), self.value(), aerosol_name, effective_mapping)


class AerosolPropertyHdf(Creator):

    pressure = param.InstanceOf(rf.Pressure)
    filename = param.Scalar(str)
    prop_name = param.Scalar(str, required=False)

    def create(self, aerosol_name=None, **kwargs):

        aerosol_file = rf.HdfFile(self.filename())

        if self.prop_name() is None:
            prop_name = aerosol_name
        else:
            prop_name = self.prop_name()

        return rf.AerosolPropertyHdf(aerosol_file, prop_name + "/Properties", self.pressure())

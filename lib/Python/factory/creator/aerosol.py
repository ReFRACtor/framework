
import os

import numpy as np

from .base import Creator, ParamPassThru
from .apriori import CreatorApriori
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

        return rf.AerosolOptical(vec_extinction, vec_properties, self.pressure(), self.relative_humidity())


class AerosolShapeGaussian(CreatorApriori):

    pressure = param.InstanceOf(rf.Pressure)
    log_space = param.Scalar(bool, default=True)

    def create(self, aerosol_name=None, **kwargs):

        if aerosol_name is None:
            raise param.ParamError("aerosol_name not supplied to creator")

        return rf.AerosolShapeGaussian(self.pressure(), self.retrieval_flag(), self.apriori(), aerosol_name, not self.log_space())

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

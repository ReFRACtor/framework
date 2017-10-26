import os

import numpy as np

from .base import Creator, ParamPassThru
from .apriori import CreatorApriori
from .. import param

from refractor import framework as rf

class AtmosphereCreator(Creator):
    "Creates an atmosphere object"

    pressure = param.InstanceOf(rf.Pressure)
    temperature = param.InstanceOf(rf.Temperature)
    altitudes = param.ObjectVector("altitude")
    absorber = param.InstanceOf(rf.Absorber)
    relative_humidity = param.InstanceOf(rf.RelativeHumidity)
    ground = param.InstanceOf(rf.Ground)
    #aerosol =
    constants = param.InstanceOf(rf.Constant)

    def create(self, **kwargs):

        pressure = self.common_store["pressure"] = self.pressure()
        temperature = self.common_store["temperature"] = self.temperature()
        altitudes = self.common_store["altitudes"] = self.altitudes()
        absorber = self.common_store["absorber"] = self.absorber()
        relative_humidity = self.relative_humidity()
        ground = self.ground()

        return rf.AtmosphereOco(absorber, pressure, temperature, relative_humidity, ground, altitudes, self.constants())

class PressureSigma(CreatorApriori):
    "Creates a PressureSigma object statisfying the AtmosphereCreator's pressure parameter"

    a_coeff = param.Array(dims=1)
    b_coeff = param.Array(dims=1)

    def create(self, **kwargs):
        # ap and flag loaded as arrays, so just get first value
        ap_psurf = self.apriori()[0]
        ret_flag = bool(self.retrieval_flag()[0])

        return rf.PressureSigma(self.a_coeff(), self.b_coeff(), ap_psurf, ret_flag)

class TemperatureMet(CreatorApriori):
    "Creates a TemperatureMet object statisfying the AtmosphereCreator's temperature parameter"
    
    met = param.InstanceOf(rf.Meteorology)
    pressure = param.InstanceOf(rf.Pressure)

    def create(self, **kwargs):
        ap_offset = self.apriori()[0]
        ret_flag = bool(self.retrieval_flag()[0])

        return rf.TemperatureMet(self.met(), self.pressure(), ap_offset, ret_flag)

class AltitudeHydrostatic(Creator):
    "Creates a AltitudeHydrostatic object statisfying the AtmosphereCreator's altitude parameter"
    
    latitude = param.ArrayWithUnit(dims=1)
    surface_height = param.ArrayWithUnit(dims=1)

    pressure = param.InstanceOf(rf.Pressure)
    temperature = param.InstanceOf(rf.Temperature)

    num_channels = param.Scalar(int)

    def create(self, **kwargs):
        # These are per channel
        latitudes = self.latitude()
        surface_heights = self.surface_height()

        altitudes = rf.vector_altitude()
        for chan_idx in range(self.num_channels()):
            chan_alt = rf.AltitudeHydrostatic(self.pressure(), self.temperature(), latitudes[chan_idx], surface_heights[chan_idx])
            altitudes.push_back(chan_alt)

        return altitudes

class ConstantForAllLevels(Creator):

    pressure = param.InstanceOf(rf.Pressure)
    value = param.Scalar()

    def create(self, **kwargs):
        return np.full(self.pressure().max_number_level, self.value())

class RelativeHumidity(Creator):

    pressure = param.InstanceOf(rf.Pressure)
    temperature = param.InstanceOf(rf.Temperature)
    absorber = param.InstanceOf(rf.Absorber)

    def create(self, **kwargs):
        return rf.RelativeHumidity(self.absorber(), self.temperature(), self.pressure())

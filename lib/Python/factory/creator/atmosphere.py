import os

import numpy as np

from .base import Creator, ParamPassThru, CreatorError
from .value import CreatorFlaggedValue
from .. import param

from refractor import framework as rf

class AtmosphereCreator(Creator):
    "Creates an atmosphere object"

    pressure = param.InstanceOf(rf.Pressure)
    temperature = param.InstanceOf(rf.Temperature)
    altitude = param.ObjectVector("altitude")
    absorber = param.InstanceOf(rf.Absorber)
    relative_humidity = param.InstanceOf(rf.RelativeHumidity)
    ground = param.InstanceOf(rf.Ground, required=False)
    aerosol = param.InstanceOf(rf.Aerosol, required=False)
    surface_temperature = param.InstanceOf(rf.SurfaceTemperature, required=False)
    constants = param.InstanceOf(rf.Constant)

    def create(self, **kwargs):

        pressure = self.common_store["pressure"] = self.pressure()
        temperature = self.common_store["temperature"] = self.temperature()
        altitude = self.common_store["altitude"] = self.altitude()
        absorber = self.common_store["absorber"] = self.absorber()
        relative_humidity = self.common_store["relative_humidity"] = self.relative_humidity()
        aerosol = self.common_store["aerosol"] = self.aerosol()
        ground = self.common_store["ground"] = self.ground()
        surf_temp = self.common_store["surface_temperature"] = self.surface_temperature()

        if surf_temp and not ground:
            raise CreatorError("Surface temperature can not be defined without ground being defined for atmosphere setup")

        if aerosol and ground and surf_temp:
            return rf.AtmosphereOco(absorber, pressure, temperature, aerosol, relative_humidity, ground, surf_temp, altitude, self.constants())
        if ground and surf_temp:
            return rf.AtmosphereOco(absorber, pressure, temperature, relative_humidity, ground, surf_temp, altitude, self.constants())
        elif aerosol and ground:
            return rf.AtmosphereOco(absorber, pressure, temperature, aerosol, relative_humidity, ground, altitude, self.constants())
        elif aerosol:
            return rf.AtmosphereOco(absorber, pressure, temperature, aerosol, relative_humidity, altitude, self.constants())
        elif ground:
            return rf.AtmosphereOco(absorber, pressure, temperature, relative_humidity, ground, altitude, self.constants())
        else:
            return rf.AtmosphereOco(absorber, pressure, temperature, relative_humidity, altitude, self.constants())


class PressureSigma(CreatorFlaggedValue):
    "Creates a PressureSigma object statisfying the AtmosphereCreator's pressure parameter"

    a_coeff = param.Array(dims=1)
    b_coeff = param.Array(dims=1)

    def create(self, **kwargs):
        # value and flag loaded as arrays, so just get first value
        psurf = self.value()[0]
        ret_flag = bool(self.retrieval_flag()[0])

        return rf.PressureSigma(self.a_coeff(), self.b_coeff(), psurf, ret_flag)

class PressureGrid(CreatorFlaggedValue):
    "Creates a PressureSigma object statisfying the AtmosphereCreator's pressure parameter"

    pressure_levels = param.Array(dims=1)

    def create(self, **kwargs):
        # ap and flag loaded as arrays, so just get first value
        psurf = self.value()[0]
        ret_flag = bool(self.retrieval_flag()[0])

        return rf.PressureSigma(self.pressure_levels(), psurf, ret_flag)

class TemperatureMet(CreatorFlaggedValue):
    "Creates a TemperatureMet object statisfying the AtmosphereCreator's temperature parameter"
    
    met = param.InstanceOf(rf.Meteorology)
    pressure = param.InstanceOf(rf.Pressure)

    def create(self, **kwargs):
        offset = self.value()[0]
        ret_flag = bool(self.retrieval_flag()[0])

        return rf.TemperatureMet(self.met(), self.pressure(), offset, ret_flag)

class SurfaceTemperature(CreatorFlaggedValue):
    "Creates a SurfaceTemperature object for use by AtmospherCreator"

    def create(self, **kwargs):

        return rf.SurfaceTemperatureDirect(self.value(), self.retrieval_flag())

class TemperatureLevelOffset(CreatorFlaggedValue):
    "Creates a TemperatureMet object statisfying the AtmosphereCreator's temperature parameter"
    
    temperature_levels = param.Array(dims=1)
    pressure = param.InstanceOf(rf.Pressure)

    def create(self, **kwargs):
        offset = self.value()[0]
        ret_flag = bool(self.retrieval_flag()[0])

        return rf.TemperatureLevelOffset(self.pressure(), self.temperature_levels(), offset, ret_flag)

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

        altitude = rf.vector_altitude()
        for chan_idx in range(self.num_channels()):
            chan_alt = rf.AltitudeHydrostatic(self.pressure(), self.temperature(), latitudes[chan_idx], surface_heights[chan_idx])
            altitude.push_back(chan_alt)

        return altitude

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

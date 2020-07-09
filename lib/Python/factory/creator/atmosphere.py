import os
import warnings

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
    rayleigh = param.InstanceOf(rf.Rayleigh, required=False)
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
        constants = self.common_store["constants"] = self.constants()
        rayleigh = self.common_store["rayleigh"] = self.rayleigh()

        if surf_temp and not ground:
            raise CreatorError("Surface temperature can not be defined without ground being defined for atmosphere setup")

        if rayleigh is None:
            warnings.warn("Rayleigh not explicitly specified using backwards compatibility component: RayleighYoung")
            rayleigh = rf.RayleighYoung(pressure, altitude, constants)

        if aerosol and ground and surf_temp:
            return rf.AtmosphereStandard(absorber, pressure, temperature, rayleigh, aerosol, relative_humidity, ground, surf_temp, altitude, constants)
        if ground and surf_temp:
            return rf.AtmosphereStandard(absorber, pressure, temperature, rayleigh, relative_humidity, ground, surf_temp, altitude, constants)
        elif aerosol and ground:
            return rf.AtmosphereStandard(absorber, pressure, temperature, rayleigh, aerosol, relative_humidity, ground, altitude, constants)
        elif aerosol:
            return rf.AtmosphereStandard(absorber, pressure, temperature, rayleigh, aerosol, relative_humidity, altitude, constants)
        elif ground:
            return rf.AtmosphereStandard(absorber, pressure, temperature, rayleigh, relative_humidity, ground, altitude, constants)
        else:
            return rf.AtmosphereStandard(absorber, pressure, temperature, rayleigh, relative_humidity, altitude, constants)


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

class TemperatureLevel(CreatorFlaggedValue):
    "Creates a TemperatureLevel object statisfying the AtmosphereCreator's temperature parameter"
    
    pressure = param.InstanceOf(rf.Pressure)

    def create(self, **kwargs):
        return rf.TemperatureLevel(self.value(), self.retrieval_flag(), self.pressure())

class TemperatureLevelOffset(CreatorFlaggedValue):
    """Creates a TemperatureLevel object statisfying the AtmosphereCreator's temperature parameter. 
    Instead of all temperature levels only an offset is retrieved."""
    
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

class AtmosphereDictCreator(Creator):
    "Creates a dictionary of atmosphere components in the correct order"

    pressure = param.InstanceOf(rf.Pressure)
    temperature = param.InstanceOf(rf.Temperature)
    absorber = param.Iterable(rf.AbsorberVmr)
    ground = param.InstanceOf(rf.Ground, required=False)
    aerosol = param.InstanceOf(rf.Aerosol, required=False)
    surface_temperature = param.InstanceOf(rf.SurfaceTemperature, required=False)

    def create(self, **kwargs):

        pressure = self.common_store["pressure"] = self.pressure()
        temperature = self.common_store["temperature"] = self.temperature()
        absorber = self.common_store["absorber"] = self.absorber()
        aerosol = self.common_store["aerosol"] = self.aerosol()
        ground = self.common_store["ground"] = self.ground()
        surf_temp = self.common_store["surface_temperature"] = self.surface_temperature()

        return  {
            'pressure': pressure,
            'temperature': temperature,
            'absorber': absorber,
            'aerosol': aerosol,
            'ground': ground,
            'surface_temperature': surf_temp,
        }

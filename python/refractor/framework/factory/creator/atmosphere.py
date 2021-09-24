import os
import warnings

from attrdict import AttrDict
import numpy as np

from .base import Creator, ParamPassThru, CreatorError
from .. import param

import refractor.framework as rf

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


class PressureSigma(Creator):
    "Creates a PressureSigma object statisfying the AtmosphereCreator's pressure parameter"

    surface_pressure = param.Scalar(float)
    a_coeff = param.Array(dims=1)
    b_coeff = param.Array(dims=1)

    def create(self, **kwargs):

        return rf.PressureSigma(self.a_coeff(), self.b_coeff(), self.surface_pressure())

class PressureGrid(Creator):
    "Creates a PressureSigma object statisfying the AtmosphereCreator's pressure parameter"

    surface_pressure = param.Scalar(float)
    pressure_levels = param.Array(dims=1)

    def create(self, **kwargs):

        return rf.PressureSigma(self.pressure_levels(), self.surface_pressure())

class SurfacePressureFromAltitude(Creator):

    pressure = param.Choice(param.InstanceOf(rf.Pressure), param.Array(dims=1))
    temperature_profile = param.Array(dims=1)
    surface_height = param.Choice(param.ArrayWithUnit(dims=1), param.DoubleWithUnit())
    latitude = param.Choice(param.ArrayWithUnit(dims=1), param.DoubleWithUnit())

    def create(self, **kwargs):
        "Calculate surface pressure from surface altitude using hydrostatic equation"

        press_in = self.pressure()

        if isinstance(press_in, rf.Pressure):
            press_obj = press_in
            press_levels = press_obj.pressure_grid().value.value
        else:
            press_levels = self.pressure()
            press_obj = rf.PressureSigma(press_levels, press_levels[-1])

        temp_levels = self.temperature_profile()
        temp_obj = rf.TemperatureLevel(temp_levels, press_obj)

        lat_in = self.latitude()
        if isinstance(lat_in, rf.DoubleWithUnit):
            lat_val = lat_in
        else:
            lat_val = lat_in[0]

        sea_level_height = rf.DoubleWithUnit(0, "m")
        alt_calc = rf.AltitudeHydrostatic(press_obj, temp_obj, lat_val, sea_level_height)

        alt_grid = np.zeros(press_levels.shape)
        for lev_idx in range(press_levels.shape[0]):
            press_val = rf.AutoDerivativeWithUnitDouble(press_obj.pressure_grid().value[lev_idx], 
                                                        press_obj.pressure_grid().units)
            alt_grid[lev_idx] = alt_calc.altitude(press_val).convert("m").value.value

        surf_height_in = self.surface_height()
        if isinstance(surf_height_in, rf.DoubleWithUnit):
            surf_height_val = self.surface_height().convert("m").value
        else:
            surf_height_val = self.surface_height()[0].convert("m").value

        # Altitude grid must be in increasing order
        surf_press = np.interp(surf_height_val, alt_grid[::-1], press_levels[::-1])

        return surf_press

class TemperatureMet(Creator):
    "Creates a TemperatureMet object statisfying the AtmosphereCreator's temperature parameter"
    
    offset = param.Scalar(float)
    met = param.InstanceOf(rf.Meteorology)
    pressure = param.InstanceOf(rf.Pressure)

    def create(self, **kwargs):

        return rf.TemperatureMet(self.met(), self.pressure(), self.offset())

class SurfaceTemperature(Creator):
    "Creates a SurfaceTemperature object for use by AtmospherCreator"

    surface_temperature = param.ArrayWithUnit(dims=1)

    def create(self, **kwargs):

        return rf.SurfaceTemperatureDirect(self.surface_temperature())

class TemperatureLevel(Creator):
    "Creates a TemperatureLevel object statisfying the AtmosphereCreator's temperature parameter"
    
    temperature_profile = param.Array(dims=1)
    pressure = param.InstanceOf(rf.Pressure)

    def create(self, **kwargs):
        return rf.TemperatureLevel(self.temperature_profile(), self.pressure())

class TemperatureLevelOffset(Creator):
    """Creates a TemperatureLevel object statisfying the AtmosphereCreator's temperature parameter. 
    Instead of all temperature levels only an offset is retrieved."""
    
    offset = param.Scalar(float)
    temperature_profile = param.Array(dims=1)
    pressure = param.InstanceOf(rf.Pressure)

    def create(self, **kwargs):

        return rf.TemperatureLevelOffset(self.pressure(), self.temperature_profile(), self.offset())

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
    altitude = param.ObjectVector("altitude")
    cached_alt_grid = None

    def create(self, **kwargs):

        pressure = self.common_store["pressure"] = self.pressure()
        temperature = self.common_store["temperature"] = self.temperature()
        absorber = self.common_store["absorber"] = self.absorber()
        aerosol = self.common_store["aerosol"] = self.aerosol()
        ground = self.common_store["ground"] = self.ground()
        surf_temp = self.common_store["surface_temperature"] = self.surface_temperature()
        altitude = self.altitude()

        def alt_grid(spec_index):
            # TODO: Check if cache invalidation necessary
            if not self.cached_alt_grid:
                # Direct translation from C++
                pres_grid = self.pressure().pressure_grid()
                # TODO: swig missing __call__/operator() for ArrayAdWithUnit to get AutoDerivativeWithUnit
                first_pres = rf.AutoDerivativeWithUnitDouble(pres_grid.value[0], pres_grid.units)
                alt_unit = self.altitude()[spec_index].altitude(first_pres).units
                alts = []
                # TODO: Why is this so slow? Actually using cached pressure?
                for (pres_ind, pres_val) in enumerate(pres_grid.value):
                    this_pres = rf.AutoDerivativeWithUnitDouble(pres_val, pres_grid.units)
                    alts.append(self.altitude()[spec_index].altitude(this_pres).value)
                alt_array_ad = rf.array_ad.np_to_array_ad(np.array(alts))
                self.cached_alt_grid = rf.ArrayAdWithUnit_double_1(alt_array_ad, alt_unit)
            return self.cached_alt_grid


        return  AttrDict({
            'pressure': pressure,
            'temperature': temperature,
            'absorber': absorber,
            'aerosol': aerosol,
            'ground': ground,
            'surface_temperature': surf_temp,
            'altitude': altitude,
            'alt_grid': alt_grid
        })

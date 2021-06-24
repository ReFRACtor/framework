from itertools import zip_longest
from collections.abc import Iterable

import numpy as np

from .base import Creator
from .. import param

import refractor.framework as rf

class SolarAbsorptionAndContinuum(Creator):

    doppler = param.Iterable(rf.SolarDopplerShift)
    absorption = param.Iterable(rf.SolarAbsorptionSpectrum)
    continuum = param.Iterable(rf.SolarContinuumSpectrum)

    def create(self, **kwargs):

        solar = []
        for doppler, absorption, continuum in zip(self.doppler(), self.absorption(), self.continuum()):
            solar.append( rf.SolarAbsorptionAndContinuum(doppler, absorption, continuum) )
        return solar

class SolarReferenceSpectrumAsciiFile(Creator):

    doppler = param.Iterable(rf.SolarDopplerShift, required=False)
    solar_data_files = param.Iterable(str)
    domain_units = param.Scalar(str)
    range_units = param.Scalar(str)
    num_channels = param.Scalar(int)

    def create(self, **kwargs):

        filenames = self.solar_data_files()
        doppler_shifts = self.doppler()

        if doppler_shifts is None:
            doppler_shifts = [None] * len(filenames)

        solar = []
        for doppler, solar_file in zip_longest(doppler_shifts, filenames, fillvalue=filenames[0]):
            file_data = np.loadtxt(solar_file)

            spec_domain = rf.SpectralDomain(file_data[:, 0], rf.Unit(self.domain_units()))
            spec_range = rf.SpectralRange(file_data[:, 1], rf.Unit(self.range_units()))
            solar_spec = rf.Spectrum(spec_domain, spec_range)
            
            ref_spec = rf.SolarReferenceSpectrum(solar_spec, doppler)

            solar.append(ref_spec)

        return solar

class SolarReferenceSpectrum(Creator):
    """Create a solar reference spectrum where all bands are specified in one file where we
    rely upon interpolation to those values to get the appropriate values"""

    doppler = param.Iterable(rf.SolarDopplerShift, required=False)

    # A single set of values for all channels or one per channel
    grid = param.Choice(param.ArrayWithUnit(dims=1), param.Iterable(rf.ArrayWithUnit_double_1))
    irradiance = param.Choice(param.ArrayWithUnit(dims=1), param.Iterable(rf.ArrayWithUnit_double_1))
    num_channels = param.Scalar(int)

    def normalize_input(self, input_type, values, num_channels):
        """Ensure we always have an array of ArrayWithUnit classes whether only one is supplied or multiple
        check the size of the iterable."""

        if isinstance(values, Iterable):
            if len(values) != num_channels:
                raise param.ParamError(f"{input_type} must have {num_channels} values if using an value per channel")
        elif isinstance(values, rf.ArrayWithUnit_double_1):
            values = [values] * num_channels
        else:
            raise param.ParamError(f"Unhandled input type for {input_type}: {values}")

        return values

    def create(self, **kwargs):

        num_channels = self.num_channels()

        doppler_shifts = self.doppler()

        if doppler_shifts is None:
            doppler_shifts = [None] * num_channels

        grid_values = self.normalize_input("grid", self.grid(), num_channels)
        irradiance_values = self.normalize_input("irradiance", self.irradiance(), num_channels)

        solar = []
        for grid, irradiance, doppler in zip(grid_values, irradiance_values, doppler_shifts):
            spec_domain = rf.SpectralDomain(grid)
            spec_range = rf.SpectralRange(irradiance)
            solar_spec = rf.Spectrum(spec_domain, spec_range)

            ref_spec = rf.SolarReferenceSpectrum(solar_spec, doppler)
            solar.append(ref_spec)

        return solar

class SolarDopplerShiftPolynomial(Creator):

    time = param.Iterable(rf.Time)
    latitude = param.ArrayWithUnit(dims=1)
    solar_zenith = param.ArrayWithUnit(dims=1)
    solar_azimuth = param.ArrayWithUnit(dims=1)
    surface_height = param.ArrayWithUnit(dims=1)
    constants = param.InstanceOf(rf.Constant)
    num_channels = param.Scalar(int)

    def create(self, **kwargs):

        do_doppler_shift = True

        doppler_shift = []
        for chan_index in range(self.num_channels()):
            doppler_shift.append( rf.SolarDopplerShiftPolynomial(self.time()[chan_index],
                                                              self.latitude()[chan_index],
                                                              self.solar_zenith()[chan_index],
                                                              self.solar_azimuth()[chan_index],
                                                              self.surface_height()[chan_index],
                                                              self.constants(),
                                                              do_doppler_shift) )
        return doppler_shift

class SolarDopplerShiftPolynomialFromL1b(Creator):

    l1b = param.InstanceOf(rf.Level1b)
    num_channels = param.Scalar(int)

    def create(self, **kwargs):

        l1b = self.l1b()

        doppler_shift = []
        for chan_idx in range(self.num_channels()):
            doppler = rf.SolarDopplerShiftPolynomial(l1b.time(chan_idx),
                                                     l1b.latitude(chan_idx),
                                                     l1b.solar_zenith(chan_idx),
                                                     l1b.solar_azimuth(chan_idx),
                                                     l1b.altitude(chan_idx))
            doppler_shift.append(doppler)

        return doppler_shift

class SolarAbsorptionTable(Creator):

    solar_data_file = param.Scalar(str)
    num_channels = param.Scalar(int)

    def create(self, **kwargs):
        solar_file = rf.HdfFile(self.solar_data_file())

        absorption = []
        for chan_index in range(self.num_channels()):
            absorption.append( rf.SolarAbsorptionTable(solar_file, "/Solar/Absorption/Absorption_%d" % (chan_index + 1)) )
        return absorption

class SolarContinuumTable(Creator):

    solar_data_file = param.Scalar(str)
    convert_from_photon = param.Scalar(bool, default=False)
    num_channels = param.Scalar(int)

    def create(self, **kwargs):
        solar_file = rf.HdfFile(self.solar_data_file())

        continuum = []
        for chan_index in range(self.num_channels()):
            continuum.append( rf.SolarContinuumTable(solar_file, "/Solar/Continuum/Continuum_%d" % (chan_index + 1), 
                self.convert_from_photon()) )
        return continuum

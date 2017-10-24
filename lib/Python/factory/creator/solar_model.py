import numpy as np

from .base import Creator
from .. import param

from refractor import framework as rf

class SolarAbsorptionAndContinuum(Creator):

    doppler = param.Iterable(rf.SolarDopplerShift)
    absorption = param.Iterable(rf.SolarAbsorptionSpectrum)
    continuum = param.Iterable(rf.SolarContinuumSpectrum)

    def create(self, **kwargs):

        solar = []
        for doppler, absorption, continuum in zip(self.doppler(), self.absorption(), self.continuum()):
            solar.append( rf.SolarAbsorptionAndContinuum(doppler, absorption, continuum) )
        return solar

class SolarDopplerShiftPolynomial(Creator):

    time = param.Iterable(rf.Time)
    latitude = param.ArrayWithUnit(dims=1)
    solar_zenith = param.ArrayWithUnit(dims=1)
    solar_azimuth = param.ArrayWithUnit(dims=1)
    altitude = param.ArrayWithUnit(dims=1)
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
                                                              self.altitude()[chan_index],
                                                              self.constants(),
                                                              do_doppler_shift) )
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

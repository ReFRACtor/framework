import math
import numpy as np

from .base import Creator
from .value import CreatorFlaggedValueMultiChannel
from .. import param
from .util import as_vector_string

from refractor import framework as rf

class AlbedoFromSignalLevel(Creator):

    polynomial_degree = param.Scalar(int, default=1)
    signal_level = param.ArrayWithUnit(dims=1)
    solar_zenith = param.ArrayWithUnit(dims=1)
    solar_strength = param.Array(dims=1)
    solar_distance = param.ArrayWithUnit(dims=1)
    stokes_coefficient = param.Array(dims=2)
    num_channels = param.Scalar(int)
    
    def create(self, **kwargs):

        signal = self.signal_level()
        sza_deg = self.solar_zenith()
        solar_strength = self.solar_strength()
        solar_distance = self.solar_distance()
        stokes_I = self.stokes_coefficient()[:, 0]

        albedo_val = np.zeros((self.num_channels(), self.polynomial_degree() + 1))

        for chan_idx in range(self.num_channels()):
            # Account for solar distance Fsun = Fsun0 / (solar_distance_meters/AU)^2
            # Create SolarDopplerShiftPolynomial so we can compute solar distance
            chan_solar_strength = solar_strength[chan_idx] / solar_distance[chan_idx].value**2
         
            # Account for stokes element for I
            chan_solar_strength = chan_solar_strength * stokes_I[chan_idx]

            sza_r = sza_deg[chan_idx].convert("rad").value
            albedo_val[chan_idx, 0] = math.pi * signal[chan_idx].value / (math.cos(sza_r) * chan_solar_strength) 

        return albedo_val

class GroundLambertian(CreatorFlaggedValueMultiChannel):

    band_reference = param.ArrayWithUnit(dims=1)
    desc_band_name = param.Iterable()

    def create(self, **kwargs):
        return rf.GroundLambertian(self.value(), self.retrieval_flag(), self.band_reference(), as_vector_string(self.desc_band_name()))

class GroundEmissivity(CreatorFlaggedValueMultiChannel):

    band_reference = param.ArrayWithUnit(dims=1)
    desc_band_name = param.Iterable()

    def create(self, **kwargs):
        return rf.GroundEmissivity(self.value(), self.retrieval_flag(), self.band_reference(), as_vector_string(self.desc_band_name()))

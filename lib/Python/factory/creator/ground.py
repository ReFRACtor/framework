import math
import numpy as np

from .base import Creator
from .value import CreatorFlaggedValue, CreatorFlaggedValueMultiChannel
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

class GroundEmissivityPolynomial(CreatorFlaggedValueMultiChannel):

    band_reference = param.ArrayWithUnit(dims=1)
    desc_band_name = param.Iterable()

    def create(self, **kwargs):
        return rf.GroundEmissivityPolynomial(self.value(), self.retrieval_flag(), self.band_reference(), as_vector_string(self.desc_band_name()))

class GroundEmissivityPiecewise(CreatorFlaggedValue):

    grid = param.ArrayWithUnit(dims=1)
    spec_win = param.InstanceOf(rf.SpectralWindowRange)

    def create(self, **kwargs):

        grid = self.grid()
        emiss_values = self.value()
        spec_win = self.spec_win()

        # Go through each grid value and only flag if the parameter falls within
        # the configured spectral ranges.

        ret_flag_compute = np.zeros(emiss_values.shape[0], dtype=bool)
        grid_conv = grid.convert_wave(spec_win.range_array.units)

        for grid_idx, grid_value in enumerate(grid_conv.value):
            for win_indexes in np.ndindex(spec_win.range_array.value.shape[:2]):
                # Make sure we compare values in order
                win_range = spec_win.range_array.value[win_indexes[0], win_indexes[1], :]
                win_beg = min(win_range)
                win_end = max(win_range)

                if grid_value >= win_beg and grid_value <= win_end:
                    ret_flag_compute[grid_idx] = True
                    break

        # Combine flags if the original retrieval flag was not all True, meaning it was the default one
        ret_flag_user = self.retrieval_flag()
        if not np.all(ret_flag_user):
            ret_flag = ret_flag_compute and ret_flag_user
        else:
            ret_flag = ret_flag_compute

        return rf.GroundEmissivityPiecewise(grid, emiss_values, ret_flag)



import math
import bisect
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

class GroundPiecewise(CreatorFlaggedValue):

    grid = param.ArrayWithUnit(dims=1)
    spec_win = param.InstanceOf(rf.SpectralWindowRange)

    def retrieval_flag(self):
        grid = self.grid()
        point_values = self.value()
        spec_win = self.spec_win()

        ret_flag_compute = np.zeros(point_values.shape[0], dtype=bool)

        for channel_idx in range(spec_win.range_array.value.shape[0]):
            for mw_idx in range(spec_win.range_array.value.shape[1]):
                win_range = spec_win.range_array[channel_idx, mw_idx, :]

                # Convert microwindows to same units as grid so we can keep the grid in sorted order
                win_val_1 = win_range[0].convert_wave(grid.units).value
                win_val_2 = win_range[1].convert_wave(grid.units).value

                if win_val_1 == win_val_2:
                    # Empty window, for instance where 0 = 0
                    continue

                # Make sure we compare values in order so the indexes are inte correct order after
                # the unit conversion above
                if win_val_1 < win_val_2:
                    index_beg = bisect.bisect_left(grid.value, win_val_1)
                    index_end = bisect.bisect_right(grid.value, win_val_2)
                else:
                    index_beg = bisect.bisect_left(grid.value, win_val_2)
                    index_end = bisect.bisect_right(grid.value, win_val_1)

                # If index beg and end equal zero then there are no grid points found for the window
                # This means the input grid is insufficient
                if(index_beg == 0 and index_end == 0):
                    raise param.ParamError("Input piecewise grid is insufficient to cover the microwindow range of {} to {} {}".format(win_val_1, win_val_2, grid.units))

                # Extend out matches one below the bisected indexes to ensure coverage
                index_beg = max(index_beg-1, 0)

                for idx in range(index_beg, index_end+1):
                    ret_flag_compute[idx] = True
 
        # Combine flags if the original retrieval flag was not all True, meaning it was the default one
        ret_flag_user = super().retrieval_flag()
        if not np.all(ret_flag_user):
            ret_flag = ret_flag_compute and ret_flag_user
        else:
            ret_flag = ret_flag_compute

        return ret_flag

class GroundEmissivityPiecewise(GroundPiecewise):

    def create(self, **kwargs):
        return rf.GroundEmissivityPiecewise(self.grid(), self.value(), self.retrieval_flag())

class GroundLambertianPiecewise(GroundPiecewise):

    def create(self, **kwargs):
        return rf.GroundLambertianPiecewise(self.grid(), self.value(), self.retrieval_flag())

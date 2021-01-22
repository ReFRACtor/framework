import math
import bisect
import numpy as np
from enum import Enum

from .base import Creator
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

class GroundLambertian(Creator):

    polynomial_coeffs = param.Array(dims=2)
    band_reference = param.ArrayWithUnit(dims=1)
    desc_band_name = param.Iterable()

    which_retrieved = param.Array(dims=2, required=False)

    def create(self, **kwargs):

        which_retrieved = self.which_retrieved()

        if which_retrieved is not None:
            mapping = rf.StateMappingAtIndexes(np.ravel(which_retrieved.astype(bool)))
        else:
            mapping = rf.StateMappingLinear()

        return rf.GroundLambertian(self.polynomial_coeffs(), self.band_reference(), as_vector_string(self.desc_band_name()), mapping)

class GroundEmissivityPolynomial(Creator):

    polynomial_coeffs = param.Array(dims=2)
    band_reference = param.ArrayWithUnit(dims=1)
    desc_band_name = param.Iterable()

    def create(self, **kwargs):
        return rf.GroundEmissivityPolynomial(self.polynomial_coeffs(), self.band_reference(), as_vector_string(self.desc_band_name()))

class GroundPiecewise(Creator):

    grid = param.ArrayWithUnit(dims=1)
    spec_win = param.InstanceOf(rf.SpectralWindowRange)

    def retrieved_indexes(self, point_values):
        grid = self.grid()
        spec_win = self.spec_win()

        ret_indexes = []

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
                    ret_indexes.append(idx)

        return ret_indexes
 
class GroundEmissivityPiecewise(GroundPiecewise):

    emissivity = param.Array(dims=1)

    def create(self, **kwargs):

        emissivity = self.emissivity()

        mapping = rf.StateMappingAtIndexes(self.retrieved_indexes(emissivity))

        return rf.GroundEmissivityPiecewise(self.grid(), emissivity, mapping)

class GroundLambertianPiecewise(GroundPiecewise):

    albedo = param.Array(dims=1)

    def create(self, **kwargs):

        albedo = self.albedo()

        mapping = rf.StateMappingAtIndexes(self.retrieved_indexes(albedo))

        return rf.GroundLambertianPiecewise(self.grid(), albedo, mapping)

class BrdfTypeOption(Enum):
    soil = 0
    vegetation = 1

class GroundBrdf(Creator):

    brdf_parameters = param.Array(dims=2)
    band_reference = param.ArrayWithUnit(dims=1)
    desc_band_name = param.Iterable()
    brdf_type = param.Choice(param.Scalar(int), param.InstanceOf(BrdfTypeOption), default=0)

    retrieve_kernel_params = param.Scalar(bool, default=False)

    def which_retrieved(self):
        ret_flag = np.ones(self.brdf_parameters().shape, dtype=bool)

        # Turn off all but weight offset and slope
        if not self.retrieve_kernel_params():
            ret_flag[:, :-2] = False

        # Flatten to shape as it would appear in the retrieval vector
        return ret_flag.ravel()

    def create(self, **kwargs):

        brdf_type = self.brdf_type()
        if isinstance(brdf_type, BrdfTypeOption):
            brdf_type = brdf_type.value

        mapping = rf.StateMappingAtIndexes(self.which_retrieved())

        if brdf_type == BrdfTypeOption.soil.value:
            return rf.GroundBrdfSoil(self.brdf_parameters(), self.band_reference(), as_vector_string(self.desc_band_name()), mapping)
        elif brdf_type == BrdfTypeOption.vegetation.value:
            return rf.GroundBrdfVeg(self.brdf_parameters(), self.band_reference(), as_vector_string(self.desc_band_name()), mapping)
        else:
            raise param.ParamError("Unknown BRDF type option: {}".format(brdf_type))

class BrdfWeightFromContinuum(Creator):

    solar_zenith = param.ArrayWithUnit(dims=1)
    observation_zenith = param.ArrayWithUnit(dims=1)
    relative_azimuth = param.ArrayWithUnit(dims=1)
    continuum_albedo = param.Array(dims=2)
    brdf_parameters = param.Array(dims=2)

    # No default to ensure this value does not get de-synced from what is used by GroundBrdf
    brdf_type = param.Choice(param.Scalar(int), param.InstanceOf(BrdfTypeOption))

    num_channels = param.Scalar(int)
 
    def create(self, **kwargs):
        sza = self.solar_zenith()
        vza = self.observation_zenith()
        azm = self.relative_azimuth()
        alb_cont = self.continuum_albedo()
        params = self.brdf_parameters()
        brdf_type = self.brdf_type()
        num_channels = self.num_channels()

        if isinstance(brdf_type, BrdfTypeOption):
            brdf_type = brdf_type.value

        if brdf_type == BrdfTypeOption.soil.value:
            brdf_class = rf.GroundBrdfSoil
        elif brdf_type == BrdfTypeOption.vegetation.value:
            brdf_class = rf.GroundBrdfVeg
        else:
            raise param.ParamError("Unknown BRDF type option: {}".format(brdf_type))

        # Remove offset and slope parameters when calling kernel_value
        for chan_idx in range(num_channels):
            kernel_params = params[chan_idx, 2:]

            alb_calc = brdf_class.kernel_value_at_params(kernel_params, sza[chan_idx].value, vza[chan_idx].value, azm[chan_idx].value)
            weight = alb_cont[chan_idx, 0] / alb_calc
            
            # Replace correct parameter index with new weight
            params[chan_idx, 5] = weight

        return params

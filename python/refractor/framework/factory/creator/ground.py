import math
import bisect
import numpy as np
from enum import Enum

from .base import Creator
from .. import param
from .util import as_vector_string

import refractor.framework as rf

class AlbedoFromSignalLevel(Creator):

    polynomial_degree = param.Scalar(int, default=1)
    signal_level = param.ArrayWithUnit(dims=1)
    solar_zenith = param.ArrayWithUnit(dims=1)
    solar_strength = param.Array(dims=1)
    solar_distance = param.Choice(param.ArrayWithUnit(dims=1), param.DoubleWithUnit())
    stokes_coefficient = param.Array(dims=2)
    num_channels = param.Scalar(int)
    
    def create(self, **kwargs):

        signal = self.signal_level()
        sza_deg = self.solar_zenith()
        solar_strength = self.solar_strength()
        solar_distance = self.solar_distance()
        stokes_I = self.stokes_coefficient()[:, 0]

        albedo_val = np.zeros((self.num_channels(), self.polynomial_degree() + 1))

        # Convert solar distance in meters to AU
        m_to_au = 149597870691

        for chan_idx in range(self.num_channels()):
            try:
                chan_solar_distance = solar_distance[chan_idx].value
            except:
                chan_solar_distance = solar_distance.value

            # Account for solar distance Fsun = Fsun0 / (solar_distance_meters/AU)^2
            # Create SolarDopplerShiftPolynomial so we can compute solar distance
            chan_solar_strength = solar_strength[chan_idx] / (chan_solar_distance / m_to_au)**2
         
            # Account for stokes element for I
            chan_solar_strength = chan_solar_strength * stokes_I[chan_idx]

            sza_r = sza_deg[chan_idx].convert("rad").value
            albedo_val[chan_idx, 0] = math.pi * signal[chan_idx].value / (math.cos(sza_r) * chan_solar_strength) 

        return albedo_val

class BrdfWeightFromSignalBase(AlbedoFromSignalLevel):

    sounding_zenith = param.ArrayWithUnit(dims=1)
    relative_azimuth = param.ArrayWithUnit(dims=1)

    brdf_parameters = param.Array(dims=2)

    def adjusted_params(self, brdf_class, **kwargs):

        alb_cont = super().create(**kwargs)

        sza_deg = self.solar_zenith()
        vza_deg = self.sounding_zenith()
        azm_deg = self.relative_azimuth()

        # Need to copy since we are editing this array in place and multiple calls
        # to this creator should return the same answer each time
        params = self.brdf_parameters().copy()

        for chan_idx in range(self.num_channels()):
            sza_chan = sza_deg[chan_idx].convert("deg").value
            vza_chan = vza_deg[chan_idx].convert("deg").value
            azm_chan = azm_deg[chan_idx].convert("deg").value

            # Extract all but the slope portion of the apriori to feed into the
            # albedo calculation function
            params_chan = params[chan_idx, :brdf_class.BRDF_WEIGHT_INTERCEPT_INDEX]

            alb_calc = brdf_class.kernel_value_at_params(params_chan, sza_chan, vza_chan, azm_chan)
            weight = alb_cont[chan_idx, 0] / alb_calc

            params[chan_idx, brdf_class.BRDF_WEIGHT_INTERCEPT_INDEX] *= weight

        return params

    def create(self, **kwargs):
        raise NotImplementedError("Use child classes BrdfWeightFromSignalVeg or BrdfWeightFromSignalSoil for this creator")

class BrdfWeightFromSignalVeg(BrdfWeightFromSignalBase):

    def create(self, **kwargs):
        return self.adjusted_params(rf.GroundBrdfVeg, **kwargs)

class BrdfWeightFromSignalSoil(BrdfWeightFromSignalBase):

    def create(self, **kwargs):
        return self.adjusted_params(rf.GroundBrdfSoil, **kwargs)

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

        # Use a set to avoid duplicate indexes
        ret_indexes = set()

        for channel_idx in range(spec_win.range_array.value.shape[0]):
            disp = spec_win.dispersion[channel_idx]

            for mw_idx in range(spec_win.range_array.value.shape[1]):
                win_range = spec_win.range_array[channel_idx, mw_idx, :]

                # Convert microwindows to same units as grid so we can keep the grid in sorted order
                if spec_win.range_array.units.name == 'sample_index':
                    # For spectral windows using sample_indexes we must look up what value is at the particular index
                    # These are one based indexes so subtract off 1
                    win_idx_1 = int(win_range[0].value)-1
                    win_idx_2 = int(win_range[1].value)-1

                    win_val_1 = rf.DoubleWithUnit(disp.pixel_grid.data[win_idx_1], disp.pixel_grid.units)
                    win_val_2 = rf.DoubleWithUnit(disp.pixel_grid.data[win_idx_2], disp.pixel_grid.units)
                else:
                    win_val_1 = win_range[0].convert_wave(grid.units).value
                    win_val_2 = win_range[1].convert_wave(grid.units).value

                # Convert to grid units
                win_val_1 = win_val_1.convert_wave(grid.units).value
                win_val_2 = win_val_2.convert_wave(grid.units).value

                if win_val_1 == win_val_2:
                    # Empty window, for instance where 0 = 0
                    continue

                # Make sure we compare values in order so the indexes are inte correct order after
                # the unit conversion above
                if win_val_1 < win_val_2:
                    index_beg = bisect.bisect_left(grid.value, win_val_1)
                    index_end = bisect.bisect_right(grid.value, win_val_2)-1
                else:
                    index_beg = bisect.bisect_left(grid.value, win_val_2)
                    index_end = bisect.bisect_right(grid.value, win_val_1)-1

                # If index beg and end equal zero then there are no grid points found for the window
                # This means the input grid is insufficient
                if(index_beg == 0 and index_end == 0):
                    raise param.ParamError("Input piecewise grid is insufficient to cover the microwindow range of {} to {} {}".format(win_val_1, win_val_2, grid.units))

                # Extend out matches one below the bisected indexes to ensure coverage
                index_beg = max(index_beg-1, 0)

                for idx in range(index_beg, index_end+1):
                    ret_indexes.add(idx)

        return list(ret_indexes)
 
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

class BrdfIndexes(object):

    def __init__(self, brdf_type):
        if isinstance(brdf_type, BrdfTypeOption):
            brdf_type = brdf_type.value

        self.brdf_type = brdf_type

        if brdf_type == BrdfTypeOption.soil.value:
            self.brdf_class = rf.GroundBrdfSoil
        elif brdf_type == BrdfTypeOption.vegetation.value:
            self.brdf_class = rf.GroundBrdfVeg
        else:
            raise param.ParamError("Unknown BRDF type option: {}".format(brdf_type))

    @property
    def kernel_indexes(self):
        kernel_indexes = (self.brdf_class.RAHMAN_KERNEL_FACTOR_INDEX,
                          self.brdf_class.RAHMAN_OVERALL_AMPLITUDE_INDEX,
                          self.brdf_class.RAHMAN_ASYMMETRY_FACTOR_INDEX,
                          self.brdf_class.RAHMAN_GEOMETRIC_FACTOR_INDEX,
                          self.brdf_class.BREON_KERNEL_FACTOR_INDEX)

        return kernel_indexes

    @property
    def weight_offset_index(self):
        return self.brdf_class.BRDF_WEIGHT_INTERCEPT_INDEX

    @property
    def weight_slope_index(self):
        return self.brdf_class.BRDF_WEIGHT_SLOPE_INDEX

    @classmethod
    def weight_indexes(cls, brdf_type):
        return (self.weight_offset_index, self.weight_slope_index)

class GroundBrdf(Creator):

    brdf_parameters = param.Array(dims=2)
    band_reference = param.ArrayWithUnit(dims=1)
    desc_band_name = param.Iterable()
    brdf_type = param.Choice(param.Scalar(int), param.InstanceOf(BrdfTypeOption), default=0)

    retrieve_kernel_params = param.Scalar(bool, default=False)

    def which_retrieved(self):
        ret_flag = np.ones(self.brdf_parameters().shape, dtype=bool)

        brdf_indexes = BrdfIndexes(self.brdf_type())

        # Turn off all but weight offset and slope
        if not self.retrieve_kernel_params():
            ret_flag[:, brdf_indexes.kernel_indexes] = False

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

        # Obtain parameter indexing consistent with BRDF class
        brdf_indexes = BrdfIndexes(self.brdf_type())

        for chan_idx in range(num_channels):
            kernel_params = params[chan_idx, brdf_indexes.kernel_indexes]

            alb_calc = brdf_class.kernel_value_at_params(kernel_params, sza[chan_idx].value, vza[chan_idx].value, azm[chan_idx].value)
            weight = alb_cont[chan_idx, 0] / alb_calc
            
            # Replace correct parameter index with new weight
            params[chan_idx, brdf_indexes.weight_offset_index] = weight

        return params

class GroundCoxmunk(Creator):

    # For coxmunk kernel
    windspeed = param.Scalar(dtype=float)
    refractive_index = param.Array(dims=1)
    windspeed_retrieved = param.Scalar(bool, default=True)

    # For optional lambertian portion
    albedo_coeffs = param.Array(dims=2, required=False)
    band_reference = param.ArrayWithUnit(dims=1, required=False)
    desc_band_name = param.Iterable(required=False)
    which_albedo_retrieved = param.Array(dims=2, required=False)

    def create(self, **kwargs):
        # Coxmunk portion
        windspeed_retrieved = self.windspeed_retrieved()

        if windspeed_retrieved is not None:
            mapping = rf.StateMappingAtIndexes(np.array([windspeed_retrieved]).astype(bool))
        else:
            mapping = rf.StateMappingLinear()

        coxmunk_obj = rf.GroundCoxmunk(self.windspeed(), self.refractive_index(), mapping)

        # Lambertian portion
        albedo_coeffs = self.albedo_coeffs()

        if albedo_coeffs is not None:
            lambertian_def = {
                'polynomial_coeffs': albedo_coeffs,
                'band_reference': self.band_reference,
                'desc_band_name': self.desc_band_name,
                'which_retrieved': self.which_albedo_retrieved,
            }

            lambertian_obj = GroundLambertian(lambertian_def).create()

            return rf.GroundCoxmunkPlusLambertian(coxmunk_obj, lambertian_obj)

        else:
            return coxmunk_obj

import numpy as np

from .base import Creator
from .. import param

from refractor import framework as rf

class ValueFromLevel1b(Creator):
    "Extracts the surface pressure value from the attached meterological file"

    l1b = param.InstanceOf(rf.Level1b)
    field = param.Scalar(str)

    list_fields = [
        'time',
    ]

    array_fields = [
        'stokes_coefficient',
    ]

    # Fields that can become ArrayWithUnit if no channel_index supplied
    array_with_unit_fields = [
        'latitude',
        'longitude',
        'sounding_zenith',
        'sounding_azimuth',
        'solar_zenith',
        'solar_azimuth',
        'relative_azimuth',
        'altitude',
        'relative_velocity',
        'signal',
        'spectral_coefficient',
    ]

    def _as_list(self, accessor):
        num_channels = self.l1b().number_spectrometer
        
        vals = []
        for chan_idx in range(num_channels):
            vals.append(accessor(chan_idx))
        return vals

    def _as_array(self, accessor):
        return np.array(self._as_list(accessor))
 
    def _as_array_with_unit(self, accessor):

        num_channels = self.l1b().number_spectrometer

        vals = []
        units = None
        for chan_idx in range(num_channels):
            chan_val = accessor(chan_idx)
            vals.append(chan_val.value)
            new_units = chan_val.units
            if units != None and new_units.name != units.name:
                raise param.ParamError("All units for L1B values must be the same to compact as an array")
            else:
                units = new_units

        val_arr = np.array(vals)

        if len(val_arr.shape) == 1:
            return rf.ArrayWithUnit_double_1(vals, units)
        if len(val_arr.shape) == 2:
            return rf.ArrayWithUnit_double_2(vals, units)
        else:
            raise param.ParamError("Unhandled return value size from L1b class")

    def create(self, channel_index=None, **kwargs):

        field_val = getattr(self.l1b(), self.field(), None)

        if field_val is None:
            return param.ParamError("Field does not exist in Level1b interface: %s" % self.field())

        if np.isscalar(field_val):
            return np.full(1, field_val)
        elif callable(field_val):
            if channel_index is not None:
                return field_val(channel_index)
            elif self.field() in self.array_with_unit_fields:
                return self._as_array_with_unit(field_val)
            elif self.field() in self.array_fields:
                return self._as_array(field_val)
            elif self.field() in self.list_fields:
                return self._as_list(field_val)
            else:
                return field_val()
        else:
            return field_val

class SolarDistanceFromL1b(Creator):

    l1b = param.InstanceOf(rf.Level1b)
    constants = param.InstanceOf(rf.Constant)
    do_doppler_shift = param.Scalar(bool, default=True)

    def create(self, **kwargs):

        l1b = self.l1b()
        num_channels = l1b.number_spectrometer

        solar_dist_vals = np.empty(num_channels)
        solar_dist_units = None
        for chan_idx in range(num_channels):
            chan_doppler_shift = \
                rf.SolarDopplerShiftPolynomial(l1b.time(chan_idx),
                                               l1b.latitude(chan_idx),
                                               l1b.solar_zenith(chan_idx),
                                               l1b.solar_azimuth(chan_idx),
                                               l1b.altitude(chan_idx),
                                               self.constants(),
                                               self.do_doppler_shift());
            chan_solar_dist = chan_doppler_shift.solar_distance
            solar_dist_vals[chan_idx] = chan_solar_dist.value
            solar_dist_units = chan_solar_dist.units

        return rf.ArrayWithUnit_double_1(solar_dist_vals, solar_dist_units)

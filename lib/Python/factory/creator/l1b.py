import numpy as np

from .base import Creator
from .. import param

from refractor import framework as rf

class ValueFromLevel1b(Creator):
    "Extracts the surface pressure value from the attached meterological file"

    l1b = param.InstanceOf(rf.Level1b)
    field = param.Scalar(str)

    # Fields that can become ArrayWithUnit if no channel_index supplied
    array_fields = [
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
    ]

    def _as_array(self, accessor):

        num_channels = self.l1b().number_spectrometer

        vals = np.empty(num_channels, dtype=float)
        units = None
        for chan_idx in range(num_channels):
            chan_val = accessor(chan_idx)
            vals[chan_idx] = chan_val.value
            new_units = chan_val.units
            if units != None and new_units.name != units.name:
                raise param.ParamError("All units for L1B values must be the same to compact as an array")
            else:
                units = new_units

        return rf.ArrayWithUnit_double_1(vals, units)

    def create(self, channel_index=None, **kwargs):

        field_val = getattr(self.l1b(), self.field())
        if np.isscalar(field_val):
            return np.full(1, field_val)
        elif callable(field_val):
            if channel_index is not None:
                return field_val(channel_index)
            elif self.field() in self.array_fields:
                return self._as_array(field_val)
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

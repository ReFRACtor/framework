import numpy as np

from .base import Creator
from .. import param

from refractor import framework as rf

class PerChannelMixin(object):

    num_channels = param.Scalar(int)

    def check_num_channels(self, check_num):
        "Issue a warning when the number of channels in some data source does not match the expected number"

        if not self.num_channels() == check_num:
            raise param.ParamError("Number of channels for data source: %d does not match expected number: %d" % (check_num, self.num_channels()))

class SpectralWindowRange(Creator):
    
    window_ranges = param.ArrayWithUnit(dims=3)
    bad_sample_mask = param.Array(dims=3, required=False)
    
    def create(self, **kwargs):

        win_ranges = self.window_ranges()
        bad_sample_mask = self.bad_sample_mask()

        if bad_sample_mask is not None:
            return rf.SpectralWindowRange(win_ranges, bad_sample_mask)
        else:
            return rf.SpectralWindowRange(win_ranges)

class UniformSpectrumSampling(Creator, PerChannelMixin):

    high_res_spacing = param.Choice(param.ArrayWithUnit(dims=1), param.DoubleWithUnit())

    def create(self, **kwargs):

        spacing_val = self.high_res_spacing()
        num_channels = self.num_channels()

        if isinstance(spacing_val, rf.DoubleWithUnit):
            # Create an ArrayWithDouble matching the number of channels used
            spacing_used = rf.ArrayWithUnit_double_1(np.full(num_channels, spacing_val.value), spacing_val.units)
        else:
            check_num_channels(spacing_val.value.shape[0])
            spacing_used = spacing_val

        return rf.SpectrumSamplingFixedSpacing(spacing_used)

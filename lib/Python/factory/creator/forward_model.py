from .base import Creator
from .. import param

from refractor import framework as rf

class SpectralWindowRange(Creator):
    
    window_ranges = param.ArrayWithUnit(dims=3)
    bad_sample_mask = param.Array(dims=3, required=False)
    
    def create(self):

        win_ranges = self.param("window_ranges")
        bad_sample_mask = self.param("bad_sample_mask")

        if bad_sample_mask is not None:
            return rf.SpectralWindowRange(win_ranges, bad_sample_mask)
        else:
            return rf.SpectralWindowRange(win_ranges)

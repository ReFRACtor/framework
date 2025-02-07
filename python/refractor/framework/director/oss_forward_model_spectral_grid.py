import numpy as np

import refractor.framework as rf

class OSSForwardModelSpectralGrid(rf.ForwardModelSpectralGrid):
    """Applies the OSS algorithm to model radiances."""

    def __init__(self, inst, swin, spectrum_sampling, oss_sd):
        print("OSS spectral grid python")
        self.inst = inst
        self.swin = swin
        self.spectrum_sampling = spectrum_sampling
        self.oss_sd = oss_sd
        super().__init__(inst, swin, spectrum_sampling)

    def high_resolution_grid(self, spec_index):
        return self.oss_sd
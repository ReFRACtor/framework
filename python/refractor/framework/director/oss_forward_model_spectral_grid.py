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

    def interpolate_spectrum(self, spec_in, spec_index):
        if spec_index < 0 or spec_index >= self.number_spectrometer:
            raise RuntimeError("Invalid spectrometer index")

        if not self.spectrum_sampling.need_interpolation(spec_index):
            return spec_in

        inst_sample_grid = self.inst.pixel_spectral_domain(spec_index)
        num_pts, = inst_sample_grid.data.shape

        spec_in_range = spec_in.spectral_range.data
        ispec_range = np.matmul(spec_in_range, self.interp)

        res_jacobs = np.empty((num_pts,
                               spec_in.spectral_range.data_ad.number_variable))

        if spec_in.spectral_range.data_ad.number_variable > 0:
            spec_in_range_jacob = spec_in.spectral_range.data_ad.jacobian
            for i in range(spec_in_range_jacob.shape[1]):
                res_jacobs[:, i] = np.matmul(spec_in_range_jacob[:, i], self.interp)
        res = rf.ArrayAd_double_1(ispec_range, res_jacobs)

        return rf.Spectrum(inst_sample_grid,
                           rf.SpectralRange(res, spec_in.spectral_range.units))
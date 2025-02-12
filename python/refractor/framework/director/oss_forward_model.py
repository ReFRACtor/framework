import numpy as np

import refractor.framework as rf


class OSSForwardModel(rf.StandardForwardModel):
    """Applies the OSS algorithm to model radiances."""

    def __init__(self):
        print("OSS python default const")
        super().__init__()

    def __init__(self, inst, swin, rt, spectrum_sampling, spectrum_effect, train_fname):
        print("OSS python full constr")
        self.train_fname = train_fname
        super().__init__(inst, swin, rt, spectrum_sampling, spectrum_effect)

    def setup_grid(self):
        print("OSS setup_grid")
        train_data = np.load(self.train_fname)
        self.rt_sd = rf.SpectralDomain(train_data['wl_hires_nm'], rf.Unit("nm"))
        self.oss_sd = rf.SpectralDomain(train_data['instr_grid_nm'], rf.Unit("nm"))
        self.oss_spectral_grid = rf.director.OSSForwardModelSpectralGrid(self.instrument,
                                                                     self.spectral_window,
                                                                     self.spectrum_sampling,
                                                                     self.oss_sd)
        self.oss_interp = train_data['OSS_mat_default']
        super().setup_grid()

    def spectral_domain(self, sensor_index):
        if not self.oss_spectral_grid:
            raise RuntimeError("setup_grid needs to be called before calling spectral_domain")
        return self.oss_spectral_grid.low_resolution_grid(sensor_index)

    def spectral_grid(self):
        return self.oss_spectral_grid

    def radiance(self, sensor_index, skip_jacobian=False):
        if not self.oss_spectral_grid:
            raise RuntimeError("setup_grid needs to be called before calling radiance")
        if not (0 <= sensor_index < self.num_channels):
            raise RuntimeError(f"Sensor index not in range [0, {self.num_channels})")

        oss_spec = self.radiative_transfer.reflectance(
            self.rt_sd,
            sensor_index,
            skip_jacobian)

        oss_spec_range = oss_spec.spectral_range.data
        highres_range = np.matmul(oss_spec_range, self.oss_interp)

        num_oss_pts, = self.oss_sd.data.shape
        oss_jacobs = np.empty((num_oss_pts,
                               oss_spec.spectral_range.data_ad.number_variable))

        if oss_spec.spectral_range.data_ad.number_variable > 0:
            oss_spec_range_jacob = oss_spec.spectral_range.data_ad.jacobian
            for i in range(oss_spec_range_jacob.shape[1]):
                # TODO: need to identify correct OSS_mat_jacrow_<sv>_linear_level_<level>
                oss_jacobs[:, i] = np.matmul(oss_spec_range_jacob[:, i], self.oss_interp)

        highres_range_ad = rf.ArrayAd_double_1(highres_range, oss_jacobs)

        highres_spec = rf.Spectrum(self.oss_sd,
                           rf.SpectralRange(highres_range_ad, oss_spec.spectral_range.units))
        self.notify_spectrum_update(highres_spec, "high_res_rt", sensor_index)

        # TODO: Confirm conversion from highres to uniform highres for spectrum correction applications
        # convolved_spec = self.apply_spectrum_corrections(highres_spec, sensor_index)

        # self.notify_spectrum_update(convolved_spec, "convolved", sensor_index)

        inst_wvl_nm = self.spectral_grid().low_resolution_grid(sensor_index).convert_wave(rf.Unit("nm"))
        inst_sample_grid = rf.SpectralDomain(inst_wvl_nm, rf.Unit("nm"))

        convolved_range = np.interp(inst_sample_grid.data,
                                    self.oss_sd.data,
                                    highres_range.data)

        num_inst_pts, = inst_sample_grid.data.shape
        inst_jacobs = np.empty((num_inst_pts,
                               highres_spec.spectral_range.data_ad.number_variable))

        if highres_spec.spectral_range.data_ad.number_variable > 0:
            highres_spec_range_jacob = highres_spec.spectral_range.data_ad.jacobian
            for i in range(highres_spec_range_jacob.shape[1]):
                inst_jacobs[:, i] = np.interp(inst_sample_grid.data,
                                    self.oss_sd.data,
                                    highres_spec_range_jacob[:, i])
        inst_range_ad = rf.ArrayAd_double_1(convolved_range, inst_jacobs)

        convolved_spec = rf.Spectrum(inst_sample_grid,
                           rf.SpectralRange(inst_range_ad, highres_spec.spectral_range.units))
        self.notify_spectrum_update(convolved_spec, "convolved", sensor_index)

        return convolved_spec

    def apply_spectrum_corrections(self, highres_spec, sensor_index):
        if not self.oss_spectral_grid:
            raise RuntimeError("setup_grid needs to be called before calling apply_spectrum_corrections")

        # highres_spec_intepolated = self.oss_spectral_grid.interpolate_spectrum(highres_spec, sensor_index)
        # self.notify_spectrum_update(highres_spec_intepolated, "high_res_interpolated", sensor_index);

        # for effect in self.spectrum_effect[sensor_index]:
        #     effect.apply_effect(highres_spec_intepolated, self.oss_spectral_grid)
        #     self.notify_spectrum_update(
        #         highres_spec_intepolated,
        #         f"high_res_spec_effect_{effect.name}",
        #         sensor_index)
        # return highres_spec_intepolated

        for effect in self.spectrum_effect[sensor_index]:
            effect.apply_effect(highres_spec, self.oss_spectral_grid)
            self.notify_spectrum_update(
                highres_spec,
                f"high_res_spec_effect_{effect.name}",
                sensor_index)
        return highres_spec

    def __str__(self):
        # TODO: If this method is used instead of print, add output from super().print()
        str_repr = (
            "OSSForwardModel:\n"
            f"  OSS training fname: {self.train_fname}\n"
            "  underlying forward model: \n"
            f"{self.print_parent()}"
        )
        return str_repr

    def desc(self):
        return self.__str__()

import numpy as np

from .base import Creator
from .value import CreatorFlaggedValue, CreatorFlaggedValueMultiChannel, CreatorValueMultiChannel
from .. import param

from refractor import framework as rf


class IlsInstrument(Creator):

    ils_half_width = param.ArrayWithUnit(dims=1)
    dispersion = param.Iterable(rf.SampleGrid)
    ils_function = param.Iterable(rf.IlsFunction)
    instrument_correction = param.ObjectVector("vector_instrument_correction")

    def create(self, **kwargs):

        # Store for use by other dependant creators
        dispersion = self.common_store["dispersion"] = self.dispersion()

        ils_vec = rf.vector_ils()
        for disp, ils_func, half_width in zip(dispersion, self.ils_function(), self.ils_half_width()):
            ils_vec.push_back(rf.IlsConvolution(disp, ils_func, half_width))
        return rf.IlsInstrument(ils_vec, self.instrument_correction())


class DispersionPolynomial(CreatorFlaggedValueMultiChannel):

    number_samples = param.Array(dims=1)
    is_one_based = param.Scalar(bool, default=False)
    num_parameters = param.Scalar(int, default=2)
    desc_band_name = param.Iterable(str)
    num_channels = param.Scalar(int)
    spec_win = param.InstanceOf(rf.SpectralWindow)

    def retrieval_flag(self):
        # Mask out non retrieved parameters
        ret_flag = super().retrieval_flag()
        ret_flag[:, slice(self.num_parameters(), ret_flag.shape[1])] = False
        return ret_flag

    def create(self, **kwargs):

        disp_coeffs = self.value()
        retrieval_flag = self.retrieval_flag()

        desc_band_name = self.desc_band_name()
        number_samples = self.number_samples()
        is_one_based = self.is_one_based()

        disp = []
        vec_disp = rf.vector_sample_grid()
        for chan_idx in range(self.num_channels()):
            chan_disp = rf.DispersionPolynomial(disp_coeffs.value[chan_idx, :], retrieval_flag[chan_idx, :], disp_coeffs.units,
                                                desc_band_name[chan_idx], int(number_samples[chan_idx]), is_one_based)
            disp.append(chan_disp)
            vec_disp.push_back(chan_disp)

        # Set dispersion into the spectral window class if it is a SpectralWindowRange
        # so it can use the dispersion to convert sample_indexes into an actual spectral range
        # This gets converted in the SpectralWindowRange spectral_bound method. The bounds
        # get used elsewhere to resolve which channels spectral points belong to. This
        # should definitely be handled in a better manner somehow since the order that Dispersion
        # gets sets into the SpectralWindowRange matters.
        spec_win = self.spec_win()
        if hasattr(spec_win, "dispersion"):
            self.spec_win().dispersion = vec_disp

        return disp


class SampleGridSpectralDomain(CreatorValueMultiChannel):

    number_samples = param.Array(dims=1)
    is_one_based = param.Scalar(bool, default=False)
    num_parameters = param.Scalar(int, default=2)
    num_channels = param.Scalar(int)
    # TODO: Need spectral_window addition same as DispersionPolynomial? (Add it + vec_sample_grid)

    def create(self, **kwargs):
        spectral_domains = self.value()

        number_samples = self.number_samples()
        is_one_based = self.is_one_based()

        sample_grid = []
        for chan_idx in range(self.num_channels()):
            chan_sample_grid = rf.SampleGridSpectralDomain(spectral_domains[chan_idx], "",
                                                           int(number_samples[chan_idx]), is_one_based)
            sample_grid.append(chan_sample_grid)

        return sample_grid


class IlsTable(Creator):

    # Supply these either as a combined 3D array or 3 different 2D arrays
    delta_lambda = param.Choice(param.Array(dims=3), param.Iterable(np.ndarray))
    response = param.Choice(param.Array(dims=3), param.Iterable(np.ndarray))
    wavenumber = param.Choice(param.Array(dims=3), param.Iterable(np.ndarray), default=None, required=False)

    dispersion = param.Iterable(rf.SampleGrid)
    hdf_band_name = param.Iterable(str)
    desc_band_name = param.Iterable(str)
    interpolate = param.Scalar(bool, default=False)
    log_space = param.Scalar(bool, default=False)
    num_channels = param.Scalar(int)

    def create(self, **kwargs):

        delta_lambda = self.delta_lambda()
        response = self.response()
        wavenumber = self.wavenumber()

        if type(delta_lambda) != type(response) or (wavenumber and type(wavenumber) != type(delta_lambda)):
            raise param.ParamError("delta_lambda, response and wavenumber (if supplied) must have the same data type")

        if isinstance(delta_lambda, np.ndarray):
            combined = True
        else:
            combined = False

        dispersion = self.dispersion()
        desc_band_name = self.desc_band_name()
        hdf_band_name = self.hdf_band_name()
        interpolate = self.interpolate()

        if self.log_space():
            IlsClass = rf.IlsTableLog
        else:
            IlsClass = rf.IlsTableLinear

        ils_func = []
        for chan_idx in range(self.num_channels()):
            if wavenumber is None:
                chan_wn = dispersion[chan_idx].pixel_grid.data
            else:
                if combined:
                    chan_wn = wavenumber[chan_idx, :]
                else:
                    chan_wn = wavenumber[chan_idx]

            if combined:
                ils_func.append(IlsClass(chan_wn, delta_lambda[chan_idx, :, :], response[chan_idx, :, :],
                                          desc_band_name[chan_idx], hdf_band_name[chan_idx], interpolate))
            else:
                ils_func.append(IlsClass(chan_wn, delta_lambda[chan_idx], response[chan_idx],
                                          desc_band_name[chan_idx], hdf_band_name[chan_idx], interpolate))

        return ils_func


class IlsGaussian(Creator):

    hwhm = param.ArrayWithUnit(dims=1)
    hdf_band_name = param.Iterable(str)
    desc_band_name = param.Iterable(str)
    num_channels = param.Scalar(int)

    def create(self, **kwargs):

        hwhm = self.hwhm()
        hdf_band_name = self.hdf_band_name()
        desc_band_name = self.desc_band_name()

        ils_func = []
        for chan_idx in range(self.num_channels()):
            ils_func.append(rf.IlsGaussian(hwhm.value[chan_idx], desc_band_name[chan_idx], hdf_band_name[chan_idx]))

        return ils_func


class InstrumentDoppler(CreatorFlaggedValue):

    num_channels = param.Scalar(int)

    def create(self, **kwargs):

        rel_vel = self.value()
        ret_flag = self.retrieval_flag()

        inst_doppler = []
        for chan_idx in range(self.num_channels()):
            inst_doppler.append(rf.InstrumentDoppler(rel_vel[chan_idx], bool(ret_flag[chan_idx])))

        return inst_doppler


class InstrumentCorrectionList(Creator):

    corrections = param.Iterable(str)
    num_channels = param.Scalar(int)

    def create(self, **kwargs):

        all_corrections = []
        for correction_name in self.corrections():
            self.register_parameter(correction_name, param.Iterable())
            all_corrections.append(self.param(correction_name))

        # Map these into an outer vector for each channel, with an inner vector for each correction
        inst_corr = rf.vector_vector_instrument_correction()
        for chan_index in range(self.num_channels()):
            per_channel_corr = rf.vector_instrument_correction()

            for correction in all_corrections:
                per_channel_corr.push_back(correction[chan_index])

            inst_corr.push_back(per_channel_corr)

        return inst_corr


class RadianceScaling(CreatorFlaggedValueMultiChannel):

    band_reference = param.ArrayWithUnit(dims=1)
    num_channels = param.Scalar(int)
    desc_band_name = param.Iterable(str)

    def create(self, **kwargs):

        params = self.value()
        band_ref = self.band_reference()
        retrieval_flag = self.retrieval_flag()
        band_name = self.desc_band_name()

        rad_scaling = []
        for chan_index in range(self.num_channels()):

            rad_scaling.append(rf.RadianceScalingSvFit(params[chan_index, :], retrieval_flag[chan_index, :], band_ref[chan_index], band_name[chan_index]))

        return rad_scaling


class ApplyInstrumentUnits(CreatorFlaggedValueMultiChannel):

    units = param.InstanceOf(rf.Unit)
    scale_factor = param.Scalar(float)
    num_channels = param.Scalar(int)

    def create(self, **kwargs):

        apply_units = []
        for chan_index in range(self.num_channels()):
            apply_units.append(rf.ApplyInstrumentUnits(self.units(), self.scale_factor()))

        return apply_units

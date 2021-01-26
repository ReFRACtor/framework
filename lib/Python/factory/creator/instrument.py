from warnings import warn

import numpy as np

from .base import Creator
from .. import param

from refractor import framework as rf


class IlsGratingInstrument(Creator):

    ils_half_width = param.ArrayWithUnit(dims=1)
    dispersion = param.Iterable(rf.SampleGrid)
    ils_function = param.Iterable(rf.IlsFunction)
    instrument_correction = param.ObjectVector("vector_instrument_correction")

    def create(self, **kwargs):

        # Store for use by other dependant creators
        dispersion = self.common_store["dispersion"] = self.dispersion()
        instrument_correction = self.instrument_correction()

        ils_vec = rf.vector_ils()
        for disp, ils_func, half_width in zip(dispersion, self.ils_function(), self.ils_half_width()):
            ils_vec.push_back(rf.IlsGrating(disp, ils_func, half_width))
        return rf.IlsInstrument(ils_vec, instrument_correction)


class DispersionPolynomial(Creator):

    polynomial_coeffs = param.ArrayWithUnit(dims=2)
    number_samples = param.Array(dims=1)
    is_one_based = param.Scalar(bool, required=False)
    spectral_variable = param.Iterable(np.ndarray, required=False)
    num_parameters = param.Scalar(int, default=2)
    desc_band_name = param.Iterable(str)
    num_channels = param.Scalar(int)
    spec_win = param.InstanceOf(rf.SpectralWindow)

    def retrieval_flag(self):
        # Mask out non retrieved parameters
        ret_flag = np.ones(self.polynomial_coeffs().value.shape, dtype=bool)
        ret_flag[:, slice(self.num_parameters(), ret_flag.shape[1])] = False
        return ret_flag

    def create(self, **kwargs):

        disp_coeffs = self.polynomial_coeffs()
        retrieval_flag = self.retrieval_flag()

        desc_band_name = self.desc_band_name()
        number_samples = self.number_samples()

        is_one_based = self.is_one_based()
        spec_var = self.spectral_variable()

        if is_one_based is not None and spec_var is not None:
            raise param.ParamError("is_one_based and spectral_variable cannot be defined at the same time")
        elif is_one_based is not None:
            warn("is_one_based is deprecated, instead use the spectral_variable parameter")

            if is_one_based:
                offset = 1
            else:
               offset = 0

            spec_var = [ np.arange(0, number_samples[c]) + offset for c in range(self.num_channels()) ]

        elif spec_var is None:
            raise param.ParamError("spectral_variable is a required parameter")

        disp = []
        vec_disp = rf.vector_sample_grid()
        for chan_idx in range(self.num_channels()):

            mapping = rf.StateMappingAtIndexes(retrieval_flag[chan_idx, :])

            chan_disp = rf.DispersionPolynomial(disp_coeffs[chan_idx, :], spec_var[chan_idx], desc_band_name[chan_idx], mapping)

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


class SampleGridSpectralDomain(Creator):

    spectral_domains = param.Iterable(rf.SpectralDomain)
    num_channels = param.Scalar(int)
    desc_band_name = param.Iterable(str)

    def create(self, **kwargs):

        spectral_domains = self.spectral_domains()
        desc_band_name = self.desc_band_name()

        sample_grid = []
        for chan_idx in range(self.num_channels()):
            chan_sample_grid = rf.SampleGridSpectralDomain(spectral_domains[chan_idx], desc_band_name[chan_idx])
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


class InstrumentDoppler(Creator):

    relative_velocity = param.ArrayWithUnit(dims=1)
    num_channels = param.Scalar(int)

    def create(self, **kwargs):

        rel_vel = self.relative_velocity()

        inst_doppler = []
        for chan_idx in range(self.num_channels()):
            inst_doppler.append(rf.InstrumentDoppler(rel_vel[chan_idx]))

        return inst_doppler


class InstrumentCorrectionList(Creator):

    corrections = param.Iterable(str)
    num_channels = param.Scalar(int)
    desc_band_name = param.Iterable(str)

    def create(self, **kwargs):
        
        channel_names = self.desc_band_name()

        # For each channel get all instrument corrections
        inst_corr = rf.vector_vector_instrument_correction()
        for chan_index in range(self.num_channels()):
            per_channel_corr = rf.vector_instrument_correction()

            for correction_name in self.corrections():
                # Dynamically request a InstrumentCorrection object for the desired name for the desired instrument channel index
                self.register_parameter(correction_name, param.InstanceOf(rf.InstrumentCorrection))
                per_channel_corr.push_back(self.param(correction_name, channel_index=chan_index, channel_name=channel_names[chan_index]))

            inst_corr.push_back(per_channel_corr)

        return inst_corr

class EmpiricalOrthogonalFunction(Creator):

    # Required
    scale_factors = param.Array(dims=1)
    eof_file = param.Scalar(str)
    hdf_group = param.Scalar(str)
    sounding_number = param.Scalar(int)
    order = param.Scalar(int)

    # Optional depending on EOF input file type
    # Uncertainty should be an ArrayAd_double_1 per channel
    scale_uncertainty = param.Scalar(bool, default=False)
    uncertainty = param.Iterable(rf.ArrayWithUnit_double_1, required=False)
    scale_to_stddev = param.Scalar(float, required=False)

    def create(self, channel_index=None, channel_name=None, **kwargs):

        scale_factors = self.scale_factors()
        sounding_number = self.sounding_number()
        order = self.order()

        scale_uncertainty = self.scale_uncertainty()
        scale_to_stddev = self.scale_to_stddev()
        uncertainty = self.uncertainty() 

        hdf_group = self.hdf_group()

        if scale_uncertainty and scale_to_stddev is None:
            raise param.ParamError("scale_to_stddev must be supplied when scale_uncertainty is True")

        if scale_uncertainty and uncertainty is None:
            raise param.ParamError("uncertainty must be supplied when scale_uncertainty is True")

        eof_hdf = rf.HdfFile(self.eof_file())
        
        if (scale_uncertainty):
            return rf.EmpiricalOrthogonalFunction(scale_factors[channel_index], eof_hdf, uncertainty[channel_index],
                    channel_index, sounding_number, order, channel_name, hdf_group, scale_to_stddev)
        else:
            eof_obj = rf.EmpiricalOrthogonalFunction(scale_factors[channel_index], 
                eof_hdf, channel_index, sounding_number, order, channel_name, hdf_group)

class RadianceScaling(Creator):

    scaling_params = param.Array(dims=2)
    band_reference = param.ArrayWithUnit(dims=1)

    def create(self, channel_index=None, channel_name=None, **kwargs):

        params = self.scaling_params()
        band_ref = self.band_reference()
        band_name = self.desc_band_name()

        return rf.RadianceScalingSvFit(params[channel_index, :], band_ref[channel_index], channel_name)


class ApplyInstrumentUnits(Creator):

    units = param.InstanceOf(rf.Unit)
    scale_factor = param.Scalar(float)

    def create(self, channel_index=None, channel_name=None, **kwargs):

       return rf.ApplyInstrumentUnits(self.units(), self.scale_factor())

import numpy as np

from .base import Creator
from .apriori import CreatorAprioriMultiChannel
from .. import param

from refractor import framework as rf

class IlsInstrument(Creator):

    ils_half_width = param.ArrayWithUnit(dims=1)
    dispersion = param.Iterable(rf.Dispersion)
    ils_function = param.Iterable(rf.IlsFunction)
    instrument_correction = param.ObjectVector("vector_instrument_correction")

    def create(self, **kwargs):

        # Store for use by other dependant creators
        dispersion = self.common_store["dispersion"] = self.dispersion()

        ils_vec = rf.vector_ils()
        for disp, ils_func, half_width in zip(dispersion, self.ils_function(), self.ils_half_width()):
            ils_vec.push_back( rf.IlsConvolution(disp, ils_func, half_width) )
        return rf.IlsInstrument(ils_vec, self.instrument_correction())

class DispersionPolynomial(CreatorAprioriMultiChannel):

    number_samples = param.Array(dims=1)
    is_one_based = param.Scalar(bool, default=False)
    num_parameters = param.Scalar(int, default=2)
    desc_band_name = param.Iterable(str)
    num_channels = param.Scalar(int)
    spec_win = param.InstanceOf(rf.SpectralWindow)

    def create(self, **kwargs):

        apriori = self.apriori()
        retrieval_flag = self.retrieval_flag()
        
        desc_band_name = self.desc_band_name()
        number_samples = self.number_samples()
        is_one_based = self.is_one_based()

        disp = []
        vec_disp = rf.vector_dispersion()
        for chan_idx in range(self.num_channels()):
            chan_disp = rf.DispersionPolynomial(apriori.value[chan_idx, :], retrieval_flag[chan_idx, :], apriori.units,
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
 
class IlsTable(Creator):

    dispersion = param.Iterable(rf.Dispersion)
    delta_lambda = param.Array(dims=3)
    response = param.Array(dims=3)
    hdf_band_name = param.Iterable(str)
    desc_band_name = param.Iterable(str)
    interpolate = param.Scalar(bool, default=False)
    log_space = param.Scalar(bool, default=False)
    num_channels = param.Scalar(int)

    def create(self, **kwargs):

        dispersion = self.dispersion()
        delta_lambda = self.delta_lambda()
        response = self.response()
        desc_band_name = self.desc_band_name()
        hdf_band_name = self.hdf_band_name()
        interpolate = self.interpolate()

        if self.log_space():
            IlsClass = rf.IlsTableLog
        else:
            IlsClass = rf.IlsTableLinear

        ils_func = []
        for chan_idx in range(self.num_channels()):

            wavenumber = dispersion[chan_idx].pixel_grid.data

            ils_func.append( IlsClass(wavenumber, delta_lambda[chan_idx, :, :], response[chan_idx, :, :], 
                                      desc_band_name[chan_idx], hdf_band_name[chan_idx], interpolate) )
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
            ils_func.append( rf.IlsGaussian(hwhm.value[chan_idx], desc_band_name[chan_idx], hdf_band_name[chan_idx]) )

        return ils_func

class InstrumentCorrectionList(Creator):

    corrections = param.Iterable(str)
    num_channels = param.Scalar(int)

    def create(self, **kwargs):

        # TBD
        # Loop over each thing described in corrections and add a vector_spectrum_effect of those items
        # Similar to how things are done by absorbers
        return rf.vector_vector_instrument_correction()

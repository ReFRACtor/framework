import numpy as np
import refractor.framework as rf


from test_support import *

# Doesn't currently work.
@skip
@require_serialize
def test_oss_forward_model_spectral_grid(sample_instrument, sample_spectral_window):

    wn_starts = [12963.0, 6182.0, 4809.0]
    wn_ends = [13171.0, 6256.0, 4882.0]
    wn_step = 0.05
    num_spec = 3
    inv_cm = rf.Unit("cm^-1")
    interpolated_spec = rf.SimpleFixedSpectrumSampling(wn_starts[0], wn_ends[0], wn_step,
                                                       wn_starts[1], wn_ends[1], wn_step,
                                                       wn_starts[2], wn_ends[2], wn_step)
    # Use non uniform sampling so that interpolated_spectrum will actually do something other than
    # return the high res spectrum
    grids = [np.linspace(wn_starts[spec_idx], wn_ends[spec_idx], 1000) for spec_idx in range(num_spec)]
    grids = [np.delete(grid, range(10, 900, 2)) for grid in grids]
    num_interp_samples, = grids[0].shape
    spec_domains = [rf.SpectralDomain(grid, inv_cm) for grid in grids]
    non_uni_spec_samp = rf.NonuniformSpectrumSampling(*spec_domains, interpolated_spec)

    # TODO: Replace with a loading of sample npz and update test data
    num_inst_grid, = sample_instrument.pixel_spectral_domain(0).data.shape
    w_vec_mult_full = np.ones((num_interp_samples, num_inst_grid))
    expected_spectral_range = 153735.

    spectral_grid = rf.director.OSSForwardModelSpectralGrid(sample_instrument,
                                                            sample_spectral_window,
                                                            non_uni_spec_samp,
                                                            w_vec_mult_full)
    for spec_idx in range(spectral_grid.number_spectrometer):
        # Interpolate from reduced high resolution grid to instrument grid
        inter_hgrid = rf.SpectralDomain(
            np.flip(spectral_grid.high_resolution_grid(spec_idx).convert_wave(inv_cm)), inv_cm)

        num_var = 3
        rad_data = np.array(range(num_interp_samples)).astype(np.double)
        rad_data_jacob = np.stack([range(num_interp_samples) for var in range(num_var)], axis=1)
        rad_data_ad = rf.ArrayAd_double_1(rad_data, rad_data_jacob)

        spec_range = rf.SpectralRange(rad_data_ad, rf.Unit("W / cm^2 / sr / cm^-1"))

        spec_in = rf.Spectrum(inter_hgrid, spec_range)

        spec_out = spectral_grid.interpolate_spectrum(spec_in, spec_idx)

        assert np.all(spec_out.spectral_range.data == expected_spectral_range)
    

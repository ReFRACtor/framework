import numpy as np
import refractor.framework as rf


from test_support import *

@require_serialize
def test_oss_forward_model_spectral_grid(sample_instrument, sample_spectral_window):
    wn_starts = [2324.7, 2324.7, 2324.7]
    wn_ends = [2335.4, 2335.4, 2335.4]

    wn_step = 0.05
    num_spec = 3
    nm = rf.Unit("nm")

    spec_sampling = rf.SimpleFixedSpectrumSampling(wn_starts[0], wn_ends[0], wn_step,
                                                   wn_starts[1], wn_ends[1], wn_step,
                                                   wn_starts[2], wn_ends[2], wn_step)
    swind = np.empty((num_spec,1,2))
    swind[:,0,0] = wn_starts
    swind[:,0,1] = wn_ends
    swind = rf.ArrayWithUnit_double_3(swind, nm)
    swin = rf.SpectralWindowRange(swind)

    oss_fname = (unit_test_data + "OSS_training.npz")
    train_data = np.load(oss_fname)
    oss_sd = rf.SpectralDomain(train_data['instr_grid_nm'], rf.Unit("nm"))
    oss_spectral_grid = rf.OSSForwardModelSpectralGrid(sample_instrument,
                                                       swin,
                                                       spec_sampling,
                                                       oss_sd)
    for spec_idx in range(oss_spectral_grid.number_spectrometer):
        high_res_dom = oss_spectral_grid.high_resolution_grid(spec_idx)
        assert np.all(high_res_dom.data == oss_sd.data)
    

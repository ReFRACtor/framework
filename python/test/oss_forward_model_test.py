from datetime import datetime

import numpy as np
import refractor.framework as rf


from test_support import *

# Doesn't currently work. There is a mismatch in the arguments sent to
# OSSForwardModel.
@skip
@require_serialize
def test_oss_forward_model(sample_instrument, sample_low_stream_rt, sample_high_stream_rt):
    wn_starts = [12950.0, 0.0, 0.0]
    wn_ends = [13190.0, 0.0, 0.0]
    wn_step = 0.05
    num_spec = 3
    inv_cm = rf.Unit("cm^-1")
    deg = rf.Unit("deg")

    swind = np.empty((3,1,2))
    swind[:,0,0] = wn_starts
    swind[:,0,1] = wn_ends
    swind = rf.ArrayWithUnit_double_3(swind, inv_cm)
    swin = rf.SpectralWindowRange(swind)

    cfg = rf.HdfFile(unit_test_data + "lua/example_static_input.h5")
    rt = rf.LsiRt(sample_low_stream_rt, sample_high_stream_rt, cfg)

    t = datetime(2006, 9 , 14, 12, 27, 22, 1000)
    rf_t = rf.Time().parse_time(t.isoformat())
    lat = rf.DoubleWithUnit(77.1828918457, deg)
    solar_zen = rf.DoubleWithUnit(74.128288269, deg)
    solar_az = rf.DoubleWithUnit(167.495071411, deg)
    elevation = rf.DoubleWithUnit(416, rf.Unit("m"))
    constant = rf.DefaultConstant()
    doppler = rf.SolarDopplerShiftPolynomial(rf_t, lat, solar_zen, solar_az, elevation, constant)

    hdf_static_input = rf.HdfFile(unit_test_data + "../../../input/common/input/l2_solar_model.h5")
    absorption = rf.SolarAbsorptionTable(hdf_static_input, "Solar/Absorption/Absorption_1")

    param = np.array([8.83596E21,                      
                      -9.48206E20,
                      -1.517E22,
                      1.74114E22,
                      -7.73485E21,
                      1.2313E21])
    param = rf.ArrayWithUnit_double_1(param, rf.Unit("ph / (s * m * m * micron)"))
    continuum = rf.SolarContinuumPolynomial(param)
    spec_effect = []
    for sidx in range(num_spec):
        spec_sp_eff = []
        solar_model = rf.SolarAbsorptionAndContinuum(doppler, absorption, continuum)
        spec_sp_eff.append(solar_model)
        spec_effect.append(spec_sp_eff)
    
    interpolated_spec = rf.SimpleFixedSpectrumSampling(wn_starts[0], wn_ends[0], wn_step,
                                                       wn_starts[1], wn_ends[1], wn_step,
                                                       wn_starts[2], wn_ends[2], wn_step);
    # Use non uniform sampling so that interpolated_spectrum will actually do something other than
    # return the high res spectrum
    grid = np.linspace(wn_starts[0], wn_ends[0], 1000)
    grid = np.delete(grid, range(10, 900, 2))
    grids = [grid, [], []]
    spec_domains = [rf.SpectralDomain(grid, inv_cm) for grid in grids]
    non_uni_spec_samp = rf.NonuniformSpectrumSampling(*spec_domains, interpolated_spec)
    
    w_vec_mult_full = np.ones((555,1016))
    fm = rf.OSSForwardModel(sample_instrument, swin, rt, non_uni_spec_samp, spec_effect, w_vec_mult_full)
    fm.setup_grid()

    # TODO: This requires setting abscodir to refractor_test_data/framework/absco
    rad = fm.radiance_all().spectral_range.data_ad
    tol = 1e-5
    assert (abs(rad.value[0] - 0.33772711) < tol) and (abs(rad.value[-1] - 2.08324092) < tol)
    assert abs(rad.jacobian[955,36] - 222.7255660975737) < tol

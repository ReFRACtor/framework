from datetime import datetime

import numpy as np
import refractor.framework as rf


from test_support import *

@require_serialize
def test_oss_forward_model(sample_low_stream_rt, sample_high_stream_rt):
    wn_starts = [2324.7, 2324.7, 2324.7]
    wn_ends = [2335.4, 2335.4, 2335.4]

    wn_step = 0.05
    num_spec = 3
    nm = rf.Unit("nm")
    deg = rf.Unit("deg")

    ils = []
    for spec_index in range(num_spec):
        sdom = rf.SpectralDomain(np.arange(wn_starts[spec_index], wn_ends[spec_index], wn_step), nm)
        sgrid = rf.SampleGridSpectralDomain(sdom, "BAND7")
        ils_identity = rf.IdentityIls(sgrid)
        ils.append(ils_identity)
    inst = rf.IlsInstrument(ils)

    swind = np.empty((3,1,2))
    swind[:,0,0] = wn_starts
    swind[:,0,1] = wn_ends
    swind = rf.ArrayWithUnit_double_3(swind, nm)
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
    
    spec_sampling = rf.SimpleFixedSpectrumSampling(wn_starts[0], wn_ends[0], wn_step,
                                                       wn_starts[1], wn_ends[1], wn_step,
                                                       wn_starts[2], wn_ends[2], wn_step);
    
    oss_fname = (unit_test_data + "OSS_training.npz")
    oss_jacobians = ["CO"] * 47

    fm = rf.OSSForwardModel(
        inst,
        swin,
        rt,
        spec_sampling,
        spec_effect,
        oss_fname,
        oss_jacobians
    )
    fm.setup_grid()

    # TODO: This requires setting abscodir to refractor_test_data/framework/absco
    rad = fm.radiance(0).spectral_range.data_ad
    tol = 1e-5
    assert (abs(rad.value[0] - 0.013275) < tol) and (abs(rad.value[-1] - 0.013277) < tol)
    assert abs(rad.jacobian[0,34] - 0.000129) < tol

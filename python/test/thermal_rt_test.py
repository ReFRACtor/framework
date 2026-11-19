# Integration tests of thermal emission through the full RT stack
# (AtmosphereStandard -> SpurrRt -> LIDORT/VLIDORT/First Order drivers),
# checking radiance Jacobians with respect to the state vector against
# finite differences. These exercise the blackbody chain rule (atmosphere
# temperature on levels, surface temperature) that driver-level tests cannot.

import numpy as np
import pytest
import refractor.framework as rf

from test_support import *

# A 2.06 micron (SC-band) wavenumber with appreciable 290 K emission
WN = 4810.0
SPEC_INDEX = 2


@pytest.fixture(scope="function")
def thermal_setup(lambertian_objects):
    '''The sample Lambertian atmosphere, rebuilt with a retrieved surface
    temperature so it can drive a thermal RT.'''
    atm = lambertian_objects.atmosphere
    sv = lambertian_objects.state_vector
    nspec = atm.number_spectrometer
    surf_temp = rf.SurfaceTemperatureDirect(
        rf.ArrayWithUnit_double_1(np.full(nspec, 290.0), rf.Unit("K")))
    sv.add_observer(surf_temp)
    x0 = np.append(np.array(sv.state), np.full(nspec, 290.0))
    sv.update_state(x0)
    tatm = rf.AtmosphereStandard(atm.absorber, atm.pressure, atm.temperature,
                                 atm.rayleigh, atm.aerosol, atm.relative_humidity,
                                 atm.ground, surf_temp, atm.altitude_obj,
                                 atm.constant)
    stokes = lambertian_objects.rt.low_stream_radiative_transfer.stokes_coefficient
    geom = [np.full(nspec, v) for v in (30.0, 0.0, 0.0)]  # sza, zen, azm
    return sv, x0, tatm, stokes, geom


def state_index(sv, prefix):
    return [i for i, n in enumerate(sv.state_vector_name) if n.startswith(prefix)][0]


def check_jacobian_fd(rt, sv, x0, indexes, rtol=5e-3):
    r = rt.stokes_and_jacobian_single_wn(WN, SPEC_INDEX, None)
    assert r.value[0] > 0
    for i in indexes:
        eps = 1e-4 * max(1.0, abs(x0[i]))
        x = x0.copy()
        x[i] += eps
        sv.update_state(x)
        rp = rt.stokes_single_wn(WN, SPEC_INDEX, None)
        sv.update_state(x0)
        fd = (rp[0] - r.value[0]) / eps
        assert r.jacobian[0, i] == pytest.approx(fd, rel=rtol), sv.state_vector_name[i]
    return r


def thermal_state_indexes(sv):
    return [state_index(sv, "H2O Scaling"),
            state_index(sv, "Surface Pressure"),
            state_index(sv, "Temperature Offset"),
            state_index(sv, "Ground Lambertian SC-Band Albedo Parm 1"),
            state_index(sv, f"Surface Temperature for sensor index {SPEC_INDEX}")]


@require_serialize
def test_lidort_thermal_jacobian(thermal_setup):
    sv, x0, tatm, stokes, (sza, zen, azm) = thermal_setup
    rt = rf.LidortRt(tatm, stokes, sza, zen, azm, False, 4, 8, False, False, True)
    check_jacobian_fd(rt, sv, x0, thermal_state_indexes(sv))


@require_serialize
def test_lidort_thermal_downwelling_jacobian(thermal_setup):
    sv, x0, tatm, stokes, (sza, zen, azm) = thermal_setup
    rt = rf.LidortRt(tatm, stokes, sza, zen, azm, False, 4, 8, False, False, True)
    up = rt.stokes_single_wn(WN, SPEC_INDEX, None)[0]
    rt.rt_driver.set_downwelling(True)
    # Surface temperature only reaches a BOA up-looking observer through
    # scattering, so it is not checked here
    r = check_jacobian_fd(rt, sv, x0, thermal_state_indexes(sv)[:4])
    assert r.value[0] != pytest.approx(up, rel=1e-3)


@require_serialize
def test_vlidort_lidort_thermal_consistency(thermal_setup):
    sv, x0, tatm, stokes, (sza, zen, azm) = thermal_setup
    lrt = rf.LidortRt(tatm, stokes, sza, zen, azm, False, 4, 8, False, False, True)
    vrt = rf.VlidortRt(tatm, stokes, sza, zen, azm, False, 4, 8, False, 1, False, True)
    idx = thermal_state_indexes(sv)
    rl = lrt.stokes_and_jacobian_single_wn(WN, SPEC_INDEX, None)
    rv = check_jacobian_fd(vrt, sv, x0, idx)
    assert rv.value[0] == pytest.approx(rl.value[0], rel=1e-6)
    np.testing.assert_allclose(rv.jacobian[0, idx], rl.jacobian[0, idx], rtol=1e-5)


@require_serialize
def test_first_order_thermal_jacobian(thermal_setup):
    sv, x0, tatm, stokes, (sza, zen, azm) = thermal_setup
    rt = rf.FirstOrderRt(tatm, stokes, sza, zen, azm, 1, 2, False, True)
    check_jacobian_fd(rt, sv, x0, thermal_state_indexes(sv))

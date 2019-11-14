#include "unit_test_support.h"
#include "atmosphere_fixture.h"
#include "optical_properties.h"

#include <blitz/array.h>

using namespace FullPhysics;
using namespace blitz;
BOOST_FIXTURE_TEST_SUITE(optical_properties, AtmosphereFixture)

BOOST_AUTO_TEST_CASE(sv_basis_jacobian)
{
    double test_wn = 13179.0;
    int test_chan = 0;

    // Check results with jacobian wrt to state vector
    OpticalProperties opt_prop_wrt_sv = OpticalProperties();
    opt_prop_wrt_sv.initialize(DoubleWithUnit(test_wn, units::inv_cm), test_chan, atm->absorber_ptr(), atm->rayleigh_ptr(), atm->aerosol_ptr());

    BOOST_CHECK_MATRIX_CLOSE_TOL(atm->optical_depth_wrt_state_vector(test_wn, test_chan).value(), opt_prop_wrt_sv.total_optical_depth().value(), 1e-10);
    BOOST_CHECK_MATRIX_CLOSE_TOL(atm->optical_depth_wrt_state_vector(test_wn, test_chan).jacobian(), opt_prop_wrt_sv.total_optical_depth().jacobian(), 1e-10);

    BOOST_CHECK_MATRIX_CLOSE_TOL(atm->single_scattering_albedo_wrt_state_vector(test_wn, test_chan).value(), opt_prop_wrt_sv.total_single_scattering_albedo().value(), 1e-10);
    BOOST_CHECK_MATRIX_CLOSE_TOL(atm->single_scattering_albedo_wrt_state_vector(test_wn, test_chan).jacobian(), opt_prop_wrt_sv.total_single_scattering_albedo().jacobian(), 1e-10);
}

BOOST_AUTO_TEST_CASE(rt_basis_jacobian)
{
    double test_wn = 13179.0;
    int test_chan = 0;

    // Check results with jacobian wrt to rt inputs
    OpticalPropertiesWrtRt opt_prop_wrt_rt = OpticalPropertiesWrtRt();
    opt_prop_wrt_rt.initialize(DoubleWithUnit(test_wn, units::inv_cm), test_chan, atm->absorber_ptr(), atm->rayleigh_ptr(), atm->aerosol_ptr());

    BOOST_CHECK_MATRIX_CLOSE_TOL(atm->optical_depth_wrt_iv(test_wn, test_chan).value(), opt_prop_wrt_rt.total_optical_depth().value(), 1e-10);
    BOOST_CHECK_MATRIX_CLOSE_TOL(atm->optical_depth_wrt_iv(test_wn, test_chan).jacobian(), opt_prop_wrt_rt.total_optical_depth().jacobian(), 1e-10);

    BOOST_CHECK_MATRIX_CLOSE_TOL(atm->single_scattering_albedo_wrt_iv(test_wn, test_chan).value(), opt_prop_wrt_rt.total_single_scattering_albedo().value(), 1e-10);
    BOOST_CHECK_MATRIX_CLOSE_TOL(atm->single_scattering_albedo_wrt_iv(test_wn, test_chan).jacobian(), opt_prop_wrt_rt.total_single_scattering_albedo().jacobian(), 1e-10);

    BOOST_CHECK_MATRIX_CLOSE_TOL(atm->intermediate_variable(test_wn, test_chan).jacobian(), opt_prop_wrt_rt.intermediate_jacobian(), 1e-10);
}

BOOST_AUTO_TEST_SUITE_END()

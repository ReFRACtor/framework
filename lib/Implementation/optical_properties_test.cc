#include "unit_test_support.h"
#include "atmosphere_fixture.h"
#include "optical_properties_wrt_input.h"
#include "optical_properties_wrt_rt.h"
#include "atmosphere_legacy.h"

#include <blitz/array.h>

using namespace FullPhysics;
using namespace blitz;
BOOST_FIXTURE_TEST_SUITE(optical_properties, AtmosphereFixture)

BOOST_AUTO_TEST_CASE(sv_basis_jacobian)
{
    double test_wn = 13179.0;
    int test_chan = 0;

    boost::shared_ptr<AtmosphereLegacy> atm_legacy
        (new AtmosphereLegacy(atm->absorber_ptr(), atm->pressure_ptr(), atm->temperature_ptr(),
                              atm->aerosol_ptr(), atm->relative_humidity_ptr(), atm->ground(), 
                              atm->altitude_ptr(), atm->constant_ptr()));

    // Check results with jacobian wrt to state vector
    OpticalPropertiesWrtInput opt_prop_wrt_sv = OpticalPropertiesWrtInput();
    opt_prop_wrt_sv.initialize(DoubleWithUnit(test_wn, units::inv_cm), test_chan, atm_legacy->absorber_ptr(), atm_legacy->rayleigh_ptr(), atm_legacy->aerosol_ptr());

    BOOST_CHECK_MATRIX_CLOSE_TOL(atm_legacy->optical_depth_wrt_state_vector(test_wn, test_chan).value(), opt_prop_wrt_sv.total_optical_depth().value(), 1e-10);
    BOOST_CHECK_MATRIX_CLOSE_TOL(atm_legacy->optical_depth_wrt_state_vector(test_wn, test_chan).jacobian(), opt_prop_wrt_sv.total_optical_depth().jacobian(), 1e-10);

    BOOST_CHECK_MATRIX_CLOSE_TOL(atm_legacy->single_scattering_albedo_wrt_state_vector(test_wn, test_chan).value(), opt_prop_wrt_sv.total_single_scattering_albedo().value(), 1e-10);
    BOOST_CHECK_MATRIX_CLOSE_TOL(atm_legacy->single_scattering_albedo_wrt_state_vector(test_wn, test_chan).jacobian(), opt_prop_wrt_sv.total_single_scattering_albedo().jacobian(), 1e-10);

    // Aerosol portion of pf
    ArrayAd<double, 2> frac_aer(opt_prop_wrt_sv.aerosol_fraction());
    ArrayAd<double, 3> aer_pf_expt = atm_legacy->aerosol_ptr()->pf_mom(test_wn, frac_aer);
    ArrayAd<double, 3> aer_pf_calc = opt_prop_wrt_sv.aerosol_phase_function_moments_portion();

    BOOST_CHECK_MATRIX_CLOSE_TOL(aer_pf_expt.value(), aer_pf_calc.value(), 1e-10);
    BOOST_CHECK_MATRIX_CLOSE_TOL(aer_pf_expt.jacobian(), aer_pf_calc.jacobian(), 1e-10);

    // Total pf which include rayleigh portion which is not exposed by Atmosphere interface
    ArrayAd<double, 3> tot_pf_expt_1 = atm_legacy->scattering_moment_wrt_state_vector(test_wn, test_chan);
    ArrayAd<double, 3> tot_pf_calc_1 = opt_prop_wrt_sv.total_phase_function_moments();

    BOOST_CHECK_MATRIX_CLOSE_TOL(tot_pf_expt_1.value(), tot_pf_calc_1.value(), 1e-10);
    BOOST_CHECK_MATRIX_CLOSE_TOL(tot_pf_expt_1.jacobian(), tot_pf_calc_1.jacobian(), 1e-10);

    // Check that phase function can be recomputed with a different number of moments and scattering
    ArrayAd<double, 3> tot_pf_expt_2 = atm_legacy->scattering_moment_wrt_state_vector(test_wn, test_chan, 200, 1);
    ArrayAd<double, 3> tot_pf_calc_2 = opt_prop_wrt_sv.total_phase_function_moments(200, 1);

    BOOST_CHECK_MATRIX_CLOSE_TOL(tot_pf_expt_2.value(), tot_pf_calc_2.value(), 1e-10);
    BOOST_CHECK_MATRIX_CLOSE_TOL(tot_pf_expt_2.jacobian(), tot_pf_calc_2.jacobian(), 1e-10);

}

BOOST_AUTO_TEST_CASE(rt_basis_jacobian)
{
    double test_wn = 13179.0;
    int test_chan = 0;

    boost::shared_ptr<AtmosphereLegacy> atm_legacy
        (new AtmosphereLegacy(atm->absorber_ptr(), atm->pressure_ptr(), atm->temperature_ptr(),
                              atm->aerosol_ptr(), atm->relative_humidity_ptr(), atm->ground(), 
                              atm->altitude_ptr(), atm->constant_ptr()));

    // Check results with jacobian wrt to rt inputs
    OpticalPropertiesWrtRt opt_prop_wrt_rt = OpticalPropertiesWrtRt();
    opt_prop_wrt_rt.initialize(DoubleWithUnit(test_wn, units::inv_cm), test_chan, atm_legacy->absorber_ptr(), atm_legacy->rayleigh_ptr(), atm_legacy->aerosol_ptr());

    BOOST_CHECK_MATRIX_CLOSE_TOL(atm_legacy->optical_depth_wrt_iv(test_wn, test_chan).value(), opt_prop_wrt_rt.total_optical_depth().value(), 1e-10);
    BOOST_CHECK_MATRIX_CLOSE_TOL(atm_legacy->optical_depth_wrt_iv(test_wn, test_chan).jacobian(), opt_prop_wrt_rt.total_optical_depth().jacobian(), 1e-10);

    BOOST_CHECK_MATRIX_CLOSE_TOL(atm_legacy->single_scattering_albedo_wrt_iv(test_wn, test_chan).value(), opt_prop_wrt_rt.total_single_scattering_albedo().value(), 1e-10);
    BOOST_CHECK_MATRIX_CLOSE_TOL(atm_legacy->single_scattering_albedo_wrt_iv(test_wn, test_chan).jacobian(), opt_prop_wrt_rt.total_single_scattering_albedo().jacobian(), 1e-10);

    BOOST_CHECK_MATRIX_CLOSE_TOL(atm_legacy->intermediate_variable(test_wn, test_chan).jacobian(), opt_prop_wrt_rt.intermediate_jacobian(), 1e-10);

    // Aerosol portion of pf
    ArrayAd<double, 2> frac_aer(opt_prop_wrt_rt.aerosol_fraction());
    ArrayAd<double, 3> aer_pf_expt = atm_legacy->aerosol_ptr()->pf_mom(test_wn, frac_aer);
    ArrayAd<double, 3> aer_pf_calc = opt_prop_wrt_rt.aerosol_phase_function_moments_portion();

    BOOST_CHECK_MATRIX_CLOSE_TOL(aer_pf_expt.value(), aer_pf_calc.value(), 1e-10);
    BOOST_CHECK_MATRIX_CLOSE_TOL(aer_pf_expt.jacobian(), aer_pf_calc.jacobian(), 1e-10);

    // Total pf which include rayleigh portion which is not exposed by Atmosphere interface
    ArrayAd<double, 3> tot_pf_expt_1 = atm_legacy->scattering_moment_wrt_iv(test_wn, test_chan);
    ArrayAd<double, 3> tot_pf_calc_1 = opt_prop_wrt_rt.total_phase_function_moments();

    BOOST_CHECK_MATRIX_CLOSE_TOL(tot_pf_expt_1.value(), tot_pf_calc_1.value(), 1e-10);
    BOOST_CHECK_MATRIX_CLOSE_TOL(tot_pf_expt_1.jacobian(), tot_pf_calc_1.jacobian(), 1e-10);

    // Check that phase function can be recomputed with a different number of moments and scattering
    ArrayAd<double, 3> tot_pf_expt_2 = atm_legacy->scattering_moment_wrt_iv(test_wn, test_chan, 200, 1);
    ArrayAd<double, 3> tot_pf_calc_2 = opt_prop_wrt_rt.total_phase_function_moments(200, 1);

    BOOST_CHECK_MATRIX_CLOSE_TOL(tot_pf_expt_2.value(), tot_pf_calc_2.value(), 1e-10);
    BOOST_CHECK_MATRIX_CLOSE_TOL(tot_pf_expt_2.jacobian(), tot_pf_calc_2.jacobian(), 1e-10);

}

BOOST_AUTO_TEST_SUITE_END()

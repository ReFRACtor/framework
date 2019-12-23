#include "unit_test_support.h"
#include "atmosphere_fixture.h"
#include "optical_properties_lsi.h"
#include "atmosphere_legacy.h"

#include <blitz/array.h>

using namespace FullPhysics;
using namespace blitz;

BOOST_FIXTURE_TEST_SUITE(optical_properties_lsi, AtmosphereFixture)

BOOST_AUTO_TEST_CASE(pack)
{
    Range ra = Range::all();

    double test_wn = 13179.0;
    int test_chan = 0;

    boost::shared_ptr<OpticalPropertiesWrtRt> opt_prop_wrt_rt(new OpticalPropertiesWrtRt());
    opt_prop_wrt_rt->initialize(DoubleWithUnit(test_wn, units::inv_cm), test_chan, atm->absorber_ptr(), atm->rayleigh_ptr(), atm->aerosol_ptr());

    ArrayAd<double, 2> packed_data = OpticalPropertiesLsi::pack(opt_prop_wrt_rt);

    BOOST_CHECK_MATRIX_CLOSE_TOL(opt_prop_wrt_rt->intermediate_jacobian(), packed_data.jacobian(), 1e-10);

    int packed_idx = 0;

    Array<double, 1> gas_od_expect_val(opt_prop_wrt_rt->gas_optical_depth_per_layer().value());
    Array<double, 1> gas_od_extract_val(packed_data.value()(ra, packed_idx));
    BOOST_CHECK_MATRIX_CLOSE_TOL(gas_od_expect_val, gas_od_extract_val, 1e-10);

    packed_idx += 1;

    Array<double, 1> ray_od_expect_val(opt_prop_wrt_rt->rayleigh_optical_depth().value());
    Array<double, 1> ray_od_extract_val(packed_data.value()(ra, packed_idx));
    BOOST_CHECK_MATRIX_CLOSE_TOL(ray_od_expect_val, ray_od_extract_val, 1e-10);

    packed_idx += 1;
    
    for(int aer_idx = 0; aer_idx < opt_prop_wrt_rt->number_aerosol_particles(); aer_idx++) {
        std::cerr << "aer_ext_idx = " << aer_idx << ", packed_idx = " << packed_idx << std::endl;
        Array<double, 1> aer_ext_expect_val(opt_prop_wrt_rt->aerosol_extinction_optical_depth_per_particle().value()(ra, aer_idx));
        Array<double, 1> aer_ext_extract_val(packed_data.value()(ra, packed_idx));
        BOOST_CHECK_MATRIX_CLOSE_TOL(aer_ext_expect_val, aer_ext_extract_val, 1e-10);

        packed_idx += 1;
    }

} 

BOOST_AUTO_TEST_CASE(unpack)
{
    Range ra = Range::all();

    double test_wn = 13179.0;
    int test_chan = 0;

    boost::shared_ptr<OpticalPropertiesWrtRt> opt_prop_wrt_rt(new OpticalPropertiesWrtRt());
    opt_prop_wrt_rt->initialize(DoubleWithUnit(test_wn, units::inv_cm), test_chan, atm->absorber_ptr(), atm->rayleigh_ptr(), atm->aerosol_ptr());

    ArrayAd<double, 2> packed_data = OpticalPropertiesLsi::pack(opt_prop_wrt_rt);

    boost::shared_ptr<AerosolOptical> aerosol = boost::dynamic_pointer_cast<AerosolOptical>(atm->aerosol_ptr());

    if(!aerosol) {
        throw Exception("Failed to convert aerosol class to AerosolOptical");
    }

    boost::shared_ptr<OpticalPropertiesLsi> opt_prop_up(new OpticalPropertiesLsi(packed_data, test_wn, aerosol, opt_prop_wrt_rt->number_gas_particles(), opt_prop_wrt_rt->number_aerosol_particles()));

    Array<double, 1> ray_od_expect_val(opt_prop_wrt_rt->rayleigh_optical_depth().value());
    Array<double, 1> ray_od_unpack_val(opt_prop_up->rayleigh_optical_depth().value());
    BOOST_CHECK_MATRIX_CLOSE_TOL(ray_od_expect_val, ray_od_unpack_val, 1e-10);

    Array<double, 2> ray_od_expect_jac(opt_prop_wrt_rt->rayleigh_optical_depth().jacobian());
    Array<double, 2> ray_od_unpack_jac(opt_prop_up->rayleigh_optical_depth().jacobian());
    BOOST_CHECK_MATRIX_CLOSE_TOL(ray_od_expect_jac, ray_od_unpack_jac, 1e-10);

    
    Array<double, 1> gas_od_expect_val(opt_prop_wrt_rt->gas_optical_depth_per_layer().value());
    Array<double, 1> gas_od_unpack_val(opt_prop_up->gas_optical_depth_per_layer().value());
    BOOST_CHECK_MATRIX_CLOSE_TOL(gas_od_expect_val, gas_od_unpack_val, 1e-10);

    Array<double, 2> gas_od_expect_jac(opt_prop_wrt_rt->gas_optical_depth_per_layer().jacobian());
    Array<double, 2> gas_od_unpack_jac(opt_prop_up->gas_optical_depth_per_layer().jacobian());
    BOOST_CHECK_MATRIX_CLOSE_TOL(gas_od_expect_jac, gas_od_unpack_jac, 1e-10);

    for(int aer_idx = 0; aer_idx < opt_prop_wrt_rt->number_aerosol_particles(); aer_idx++) {
        Array<double, 1> aer_ext_expect_val(opt_prop_wrt_rt->aerosol_extinction_optical_depth_per_particle().value()(ra, aer_idx));
        Array<double, 1> aer_ext_unpack_val(opt_prop_up->aerosol_extinction_optical_depth_per_particle().value()(ra, aer_idx));
        BOOST_CHECK_MATRIX_CLOSE_TOL(aer_ext_expect_val, aer_ext_unpack_val, 1e-10);

        Array<double, 2> aer_ext_expect_jac(opt_prop_wrt_rt->aerosol_extinction_optical_depth_per_particle().jacobian()(ra, aer_idx, ra));
        Array<double, 2> aer_ext_unpack_jac(opt_prop_up->aerosol_extinction_optical_depth_per_particle().jacobian()(ra, aer_idx, ra));
        BOOST_CHECK_MATRIX_CLOSE_TOL(aer_ext_expect_jac, aer_ext_unpack_jac, 1e-10);
    }

    for(int aer_idx = 0; aer_idx < opt_prop_wrt_rt->number_aerosol_particles(); aer_idx++) {
        Array<double, 1> aer_sca_expect_val(opt_prop_wrt_rt->aerosol_scattering_optical_depth_per_particle().value()(ra, aer_idx));
        Array<double, 1> aer_sca_unpack_val(opt_prop_up->aerosol_scattering_optical_depth_per_particle().value()(ra, aer_idx));
        BOOST_CHECK_MATRIX_CLOSE_TOL(aer_sca_expect_val, aer_sca_unpack_val, 1e-10);

        Array<double, 2> aer_sca_expect_jac(opt_prop_wrt_rt->aerosol_scattering_optical_depth_per_particle().jacobian()(ra, aer_idx, ra));
        Array<double, 2> aer_sca_unpack_jac(opt_prop_up->aerosol_scattering_optical_depth_per_particle().jacobian()(ra, aer_idx, ra));
        BOOST_CHECK_MATRIX_CLOSE_TOL(aer_sca_expect_jac, aer_sca_unpack_jac, 1e-10);
    }

} 

BOOST_AUTO_TEST_CASE(intermediate_variable)
{
    double test_wn = 13179.0;
    int test_chan = 0;

    boost::shared_ptr<AtmosphereLegacy> atm_legacy
        (new AtmosphereLegacy(atm->absorber_ptr(), atm->pressure_ptr(), atm->temperature_ptr(),
                              atm->aerosol_ptr(), atm->relative_humidity_ptr(), atm->ground(), 
                              atm->altitude_ptr(), atm->constant_ptr()));

    ArrayAd<double, 2> intermediate_v(atm_legacy->intermediate_variable(test_wn, test_chan));

    boost::shared_ptr<OpticalPropertiesWrtRt> opt_prop_wrt_rt(new OpticalPropertiesWrtRt());
    opt_prop_wrt_rt->initialize(DoubleWithUnit(test_wn, units::inv_cm), test_chan, atm->absorber_ptr(), atm->rayleigh_ptr(), atm->aerosol_ptr());

    ArrayAd<double, 2> packed_data = OpticalPropertiesLsi::pack(opt_prop_wrt_rt);

    BOOST_CHECK_MATRIX_CLOSE_TOL(intermediate_v.value(), packed_data.value(), 1e-10);
    BOOST_CHECK_MATRIX_CLOSE_TOL(intermediate_v.jacobian(), packed_data.jacobian(), 1e-10);

}

BOOST_AUTO_TEST_SUITE_END()

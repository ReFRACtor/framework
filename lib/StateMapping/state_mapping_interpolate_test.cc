#include "state_mapping_interpolate.h"

#include "unit_test_support.h"
#include "state_mapping_fixture.h"

using namespace FullPhysics;
using namespace blitz;

BOOST_FIXTURE_TEST_SUITE(state_mapping_interpolate, StateMappingInterpFixture)

BOOST_AUTO_TEST_CASE(linear_linear)
{
    StateMappingInterpolateLinearLinear map_interp(press_to, press_from);
    ArrayAd<double, 1> vals_to = map_interp.fm_view(vals_from_ad);

    Array<double, 1> expt_vals(10);
    expt_vals = 1.0,  2.0,  3.0,  4.0,  5.0,  6.0,  7.0,  8.0,  9.0, 10.0;

    BOOST_CHECK_MATRIX_CLOSE_TOL(vals_to.value(), expt_vals, 1e-7);
}

BOOST_AUTO_TEST_CASE(linear_log)
{
    StateMappingInterpolateLinearLog map_interp(press_to, press_from);
    ArrayAd<double, 1> vals_to = map_interp.fm_view(vals_from_ad);

    Array<double, 1> expt_vals(10);
    expt_vals = 1.0       ,  1.68851031,  2.85106706,  3.87296112,  4.893161,
                5.93554016,  6.91282185,  7.97262823,  8.92895752, 10.0;

    BOOST_CHECK_MATRIX_CLOSE_TOL(vals_to.value(), expt_vals, 1e-7);
}

BOOST_AUTO_TEST_CASE(log_log)
{
    StateMappingInterpolateLogLog map_interp(press_to, press_from);
    ArrayAd<double, 1> vals_to = map_interp.fm_view(vals_from_ad);

    Array<double, 1> expt_vals(10);
    expt_vals = 1.0,  2.0,  3.0,  4.0,  5.0,  6.0,  7.0,  8.0,  9.0, 10.0;

    BOOST_CHECK_MATRIX_CLOSE_TOL(vals_to.value(), expt_vals, 1e-7);
}

BOOST_AUTO_TEST_CASE(log_linear)
{
    StateMappingInterpolateLogLinear map_interp(press_to, press_from);
    ArrayAd<double, 1> vals_to = map_interp.fm_view(vals_from_ad);

    Array<double, 1> expt_vals(10);
    expt_vals = 1.0       ,  2.32318716,  3.09720203,  4.13803403,  5.09237651,
                6.07086629,  7.0822217 ,  8.03025399,  9.0699554 , 10.0;

    BOOST_CHECK_MATRIX_CLOSE_TOL(vals_to.value(), expt_vals, 1e-7);
}

BOOST_AUTO_TEST_CASE(serialize)
{
  if(!have_serialize_supported())
    return;

    boost::shared_ptr<StateMappingInterpolateLinearLinear> map_interp_orig =
        boost::make_shared<StateMappingInterpolateLinearLinear>(press_to, press_from);

    std::string serial_str = serialize_write_string(map_interp_orig);

    boost::shared_ptr<StateMappingInterpolateLinearLinear> map_interp_restore =
        serialize_read_string<StateMappingInterpolateLinearLinear>(serial_str);
 
    ArrayAd<double, 1> vals_expt = map_interp_orig->fm_view(vals_from_ad);
    ArrayAd<double, 1> vals_restore = map_interp_restore->fm_view(vals_from_ad);

    BOOST_CHECK_MATRIX_CLOSE_TOL(vals_restore.value(), vals_expt.value(), 1e-7);
 
}

BOOST_AUTO_TEST_SUITE_END()

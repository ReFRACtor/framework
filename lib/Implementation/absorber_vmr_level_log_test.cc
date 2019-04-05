#include "absorber_vmr_level_log.h"
#include "met_data_fixture.h"
#include "unit_test_support.h"

using namespace FullPhysics;
using namespace blitz;

BOOST_FIXTURE_TEST_SUITE(absorber_vmr_level_log, MetDataFixture)

BOOST_AUTO_TEST_CASE(basic)
{
    Array<bool, 1> flag(20);
    flag = true;
    Array<double, 1> vmr(20);
    vmr = 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20;
    AbsorberVmrLevelLog avmr(pressure, vmr, flag, "CO2");

    for(int i = 0; i < pressure->pressure_grid().rows(); ++i) {
        BOOST_CHECK_CLOSE(avmr.volume_mixing_ratio(pressure->pressure_grid()(i).value).value(), vmr(i), 1e-8);
    }

    for(int j = 0; j < vmr.rows(); j++) {
        BOOST_CHECK_CLOSE(avmr.coefficient_unmapped().value()(j), log(vmr(j)), 1e-8);
        BOOST_CHECK_CLOSE(log(avmr.coefficient().value()(j)), log(vmr(j)), 1e-8);
    }
}

BOOST_AUTO_TEST_SUITE_END()

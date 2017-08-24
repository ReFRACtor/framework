#include "absorber_vmr_fixed_level.h"
#include "met_data_fixture.h"
#include "unit_test_support.h"

using namespace FullPhysics;
using namespace blitz;

BOOST_FIXTURE_TEST_SUITE(absorber_vmr_fixed_level, MetDataFixture)

BOOST_AUTO_TEST_CASE(basic)
{
    boost::shared_ptr<PressureLevelInput> press_level(new PressureLevelInput(met_data->pressure_levels()));
    Array<bool, 1> flag(20);
    flag = true;
    Array<double, 1> vmr(20);
    vmr = 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20;
    AbsorberVmrFixedLevel avmr(pressure, press_level, flag, vmr, "CO2");
    BOOST_CHECK_MATRIX_CLOSE(avmr.volume_mixing_ratio_level(), vmr);
    Array<double, 1> grid_expect(19);
    grid_expect = 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,43.323824860025439;
    for(int i = 0; i < grid_expect.rows(); ++i)
        BOOST_CHECK_CLOSE(avmr.volume_mixing_ratio(pressure->pressure_grid()(i).value).value(), grid_expect(i), 1e-8);
}

BOOST_AUTO_TEST_SUITE_END()

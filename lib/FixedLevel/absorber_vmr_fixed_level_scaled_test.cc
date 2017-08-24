#include "absorber_vmr_fixed_level_scaled.h"
#include "met_data_fixture.h"
#include "unit_test_support.h"

using namespace FullPhysics;
using namespace blitz;

BOOST_FIXTURE_TEST_SUITE(absorber_vmr_fixed_level_scaled, MetDataFixture)

BOOST_AUTO_TEST_CASE(basic)
{
    boost::shared_ptr<PressureLevelInput> pressure_levels(new PressureLevelInput(met_data->pressure_levels()));
    Array<double, 1> vmr(20);
    vmr = 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20;
    double scale = 1.0;
    AbsorberVmrFixedLevelScaled avmr(pressure,
                                     pressure_levels,
                                     vmr, true, scale, "CO2");
    Array<double, 1> grid_expect(19);
    grid_expect = 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,43.323824860025439;
    for(int i = 0; i < grid_expect.rows(); ++i) {
        BOOST_CHECK_CLOSE(avmr.volume_mixing_ratio(pressure->pressure_grid()(i).value).value(), grid_expect(i), 1e-8);
    }

}

BOOST_AUTO_TEST_SUITE_END()

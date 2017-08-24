#include "unit_test_support.h"
#include "pressure_fixed_level_output.h"
#include "output_hdf.h"
#include "met_data_fixture.h"

using namespace FullPhysics;
using namespace blitz;
BOOST_FIXTURE_TEST_SUITE(pressure_fixed_level_output, MetDataFixture)

BOOST_AUTO_TEST_CASE(basic)
{
    Array<double, 1> met_press = met_data->pressure_levels();
    boost::shared_ptr<PressureLevelInput> pinput(new PressureLevelInput(met_press));
    boost::shared_ptr<PressureFixedLevel> p_fixed(new PressureFixedLevel(false, pinput, met_press(met_press.rows()-1)));

    boost::shared_ptr<StateVector> sv(new StateVector);
    PressureFixedLevelOutput po(p_fixed, sv);
    boost::shared_ptr<OutputHdf> out(new OutputHdf("pressure_output.h5", 20, 112, 5, 3));
    add_file_to_cleanup("pressure_output.h5");
    po.register_output(out);

    // Simple test, we just make sure that we can write output. All the
    // actual value calculation is checked in pressure unit test.

    out->write();
    BOOST_CHECK(true);
}

BOOST_AUTO_TEST_SUITE_END()



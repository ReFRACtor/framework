#include "pressure_fixed_level.h"
#include "global_fixture.h"
#include "unit_test_support.h"

using namespace FullPhysics;
using namespace blitz;

BOOST_FIXTURE_TEST_SUITE(pressure_fixed_level, GlobalFixture)

BOOST_AUTO_TEST_CASE(basic)
{
    StateVector sv;
    Array<double, 1> press_levels(3);
    press_levels = 1, 2, 3;
    double psurf = 2.5;
    boost::shared_ptr<PressureLevelInput> pinput(new PressureLevelInput(press_levels));
    Array<double, 1> press_grid_expect(3);
    press_grid_expect = 1, 2, 2.5;
    PressureFixedLevel pnone(false, pinput, psurf);
    BOOST_CHECK_CLOSE(pnone.surface_pressure().value.value(), psurf, 1e-4);
    BOOST_CHECK_MATRIX_CLOSE(pnone.pressure_grid().value.value(), 
                             press_grid_expect);
    sv.add_observer(pnone);
    Array<double, 1> x(1);
    x = 1.5;
    sv.update_state(x);
    BOOST_CHECK_CLOSE(pnone.surface_pressure().value.value(), psurf, 1e-4);
    BOOST_CHECK_MATRIX_CLOSE(pnone.pressure_grid().value.value(), 
                             press_grid_expect);
    sv.remove_observer(pnone);
    PressureFixedLevel pone(true, pinput, psurf);
    BOOST_CHECK_CLOSE(pone.surface_pressure().value.value(), psurf, 1e-4);
    BOOST_CHECK_MATRIX_CLOSE(pone.pressure_grid().value.value(), 
                             press_grid_expect);
    sv.add_observer(pone);
    sv.update_state(x);
    press_grid_expect.resize(2);
    press_grid_expect = 1, 1.5;
    BOOST_CHECK_CLOSE(pone.surface_pressure().value.value(), 1.5, 1e-4);
    BOOST_CHECK_MATRIX_CLOSE(pone.pressure_grid().value.value(), 
                             press_grid_expect);
}

BOOST_AUTO_TEST_SUITE_END()

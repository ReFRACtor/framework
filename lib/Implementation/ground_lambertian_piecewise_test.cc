#include "ground_lambertian_piecewise.h"
#include "unit_test_support.h"
#include "state_vector.h"

using namespace FullPhysics;
using namespace blitz;

BOOST_FIXTURE_TEST_SUITE(ground_lambertian_piecewise, GlobalFixture)

BOOST_AUTO_TEST_CASE(basic)
{
    IfstreamCs lamb_input(test_data_dir() + "in/lambertian/piecewise_lambertian");
    Array<double, 2> lamb_data;
    lamb_input >> lamb_data;

    ArrayWithUnit<double, 1> grid(lamb_data(Range::all(), 0), units::inv_cm);
    Array<double, 1> values(lamb_data(Range::all(), 1));
    Array<bool, 1> flag(values.shape());
    flag = true;

    auto emiss = GroundLambertianPiecewise(grid, values, flag);

    // Test mid point between each pair of wavenumbers
    for (int tst_idx = 0; tst_idx < values.rows() - 1; tst_idx++) {
        double wn1 = lamb_data(tst_idx, 0);
        double wn2 = lamb_data(tst_idx+1, 0);

        double tst_wn = (wn1 + wn2) / 2;
        double expt_value = (values(tst_idx) * (wn2 - tst_wn) + values(tst_idx+1) * (tst_wn - wn1)) / (wn2 - wn1);

        BOOST_CHECK_CLOSE(expt_value, emiss.value_at_point(DoubleWithUnit(tst_wn, units::inv_cm)).value(), 1e-10);
    }
}

BOOST_AUTO_TEST_SUITE_END()

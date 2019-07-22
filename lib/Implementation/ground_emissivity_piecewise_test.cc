#include "ground_emissivity_piecewise.h"
#include "unit_test_support.h"
#include "state_vector.h"

using namespace FullPhysics;
using namespace blitz;

BOOST_FIXTURE_TEST_SUITE(ground_emissivity_piecewise, GlobalFixture)

BOOST_AUTO_TEST_CASE(basic)
{
    IfstreamCs emiss_input(test_data_dir() + "in/emissivity/piecewise_emissivity");
    Array<double, 2> emiss_data;
    emiss_input >> emiss_data;

    ArrayWithUnit<double, 1> grid(emiss_data(Range::all(), 0), units::inv_cm);
    Array<double, 1> values(emiss_data(Range::all(), 1));
    Array<bool, 1> flag(values.shape());
    flag = true;

    auto emiss = GroundEmissivityPiecewise(grid, values, flag);

    // Test mid point between each pair of wavenumbers
    for (int tst_idx = 0; tst_idx < values.rows() - 1; tst_idx++) {
        double wn1 = emiss_data(tst_idx, 0);
        double wn2 = emiss_data(tst_idx+1, 0);

        double tst_wn = (wn1 + wn2) / 2;
        double expt_value = (values(tst_idx) * (wn2 - tst_wn) + values(tst_idx+1) * (tst_wn - wn1)) / (wn2 - wn1);

        BOOST_CHECK_CLOSE(expt_value, emiss.value_at_wavenumber(tst_wn).value(), 1e-10);
    }
}

BOOST_AUTO_TEST_SUITE_END()

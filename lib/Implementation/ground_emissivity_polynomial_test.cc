#include "ground_emissivity_polynomial.h"
#include "unit_test_support.h"
#include "state_vector.h"

using namespace FullPhysics;
using namespace blitz;
BOOST_FIXTURE_TEST_SUITE(ground_emissivity_polynomial, GlobalFixture)

BOOST_AUTO_TEST_CASE(basic)
{
    int num_channels = 1;
    Array<double, 2> spec_coeffs(num_channels, 2);
    spec_coeffs = 0.95, 1e-4;
    
    Array<bool, 2> flag(spec_coeffs.shape());
    flag = true;

    Array<double, 1> ref_point_values(num_channels);
    ref_point_values = 900;
    ArrayWithUnit<double, 1> ref_points(ref_point_values, Unit("cm^-1"));

    std::vector<std::string> band_names = { "Channel_1" };

    auto emiss = GroundEmissivityPolynomial(spec_coeffs, flag, ref_points, band_names);

    BOOST_CHECK_CLOSE(emiss.emissivity(DoubleWithUnit(655, units::inv_cm), 0).value(), 0.9255, 1e-8);
    BOOST_CHECK_CLOSE(emiss.emissivity(DoubleWithUnit(900, units::inv_cm), 0).value(), 0.95, 1e-8);
    BOOST_CHECK_CLOSE(emiss.emissivity(DoubleWithUnit(1090, units::inv_cm), 0).value(), 0.969, 1e-8);

    // Check jacobians, need to attach to a StateVector
    StateVector sv;
    sv.add_observer(emiss);
    sv.update_state(spec_coeffs(0, Range::all()));

    ArrayAd<double, 1> surf_param = emiss.surface_parameter(1090, 0);

    BOOST_CHECK_EQUAL(surf_param.value().rows(), 1);
  
    BOOST_CHECK_CLOSE(surf_param.value()(0), 0.969, 1e-6);
  
    // Slope jacobian should just be input_wn - ref_wn
    BOOST_CHECK_CLOSE(surf_param.jacobian()(0, 0), 1, 1e-8);
    BOOST_CHECK_CLOSE(surf_param.jacobian()(0, 1), 1090-900, 1e-6);
}

BOOST_AUTO_TEST_SUITE_END()




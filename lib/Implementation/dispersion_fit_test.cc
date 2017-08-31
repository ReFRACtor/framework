#include "unit_test_support.h"
#include "configuration_fixture.h"
#include "dispersion_fit.h"

using namespace FullPhysics;
using namespace blitz;

BOOST_FIXTURE_TEST_SUITE(dispersion_fit, ConfigurationFixture)

BOOST_AUTO_TEST_CASE(basic)
{
  DoubleWithUnit solar_line_location(0.77010968653012513, units::micron);
  DoubleWithUnit solar_line_width(pow(1.6e-5,2), units::micron);
  DoubleWithUnit search_width(1.682e-4, units::micron);
  DoubleWithUnit ils_offset(0, units::micron);

  ArrayWithUnit<double, 1> offset_scaling; 
  offset_scaling.value.resize(3);
  offset_scaling.value = 1, 0.48, 0.38;
  offset_scaling.units = Unit("cm^-1");

  Array<double, 2> disp_expt(3,5);
  disp_expt = 
    7.5740854222e-01, 1.7355000000e-05, -2.7849100000e-09, 0.0000000000e+00, 0.0000000000e+00, 
    1.5899286603e+00, 3.6268600000e-05, -5.7611500000e-09, -4.3324500000e-14, 0.0000000000e+00, 
    2.0412031060e+00, 4.7006100000e-05, -8.3350800000e-09, 8.1447200000e-13, -1.7970000000e-16;

  Array<double, 2> disp_initial(disp_expt.shape());
  disp_initial = disp_expt;
  disp_initial(Range::all(), 0) += 1e-4 * offset_scaling.value(Range::all());;

  DispersionFit disp_fit = DispersionFit(config_level_1b);
  Array<double, 2> disp_sol = disp_fit.fit(disp_initial, solar_line_location, solar_line_width, search_width, ils_offset, offset_scaling);

  BOOST_CHECK_MATRIX_CLOSE_TOL(disp_expt(Range::all(), 0), disp_sol(Range::all(), 0), 2e-6);

  if(false) {
      Array<double, 1> diff(disp_sol(Range::all(), 0)-disp_expt(Range::all(), 0)) ;
      std::cout << std::setprecision(10) << std::scientific << "Dispersion solution: " << std::endl << disp_sol << std::endl
        << "Diff: " << diff << std::endl;
  }

}

BOOST_AUTO_TEST_SUITE_END()


#include "aerosol_property_rh_hdf.h"
#include "atmosphere_fixture.h"
#include "pressure_sigma.h"
#include "unit_test_support.h"

using namespace FullPhysics;
using namespace blitz;

BOOST_FIXTURE_TEST_SUITE(aerosol_property_rh_hdf, AtmosphereFixture)

BOOST_AUTO_TEST_CASE(basic)
{
  HdfFile h(test_data_dir() + "in/aerosol/l2_merra_aerosol_RH.h5");
  AerosolPropertyRhHdf a(h, "SS/Properties", 
                         atm->pressure_ptr(), atm->relative_humidity_ptr());
  BOOST_CHECK_CLOSE(a.extinction_coefficient_each_layer(13000).value()(0), 
                    96.103672325921167, 1e-8);
  BOOST_CHECK_CLOSE(a.extinction_coefficient_each_layer(13000).value()(17), 
                    262.38294394587894, 1e-8);
  BOOST_CHECK_CLOSE(a.scattering_coefficient_each_layer(13000).value()(0), 
                    96.092536042985984, 1e-8);
  BOOST_CHECK_CLOSE(a.scattering_coefficient_each_layer(13000).value()(17), 
                    262.38294394587894, 1e-8);
  // Phase function is large, so we just check the first couple of moments.
  Array<double, 2> pf_expect(2, 6);
  pf_expect =
    1, 0.005262347806, 0.9134154949, 0.01427147261, 1, 0.9134154949,
    2.379920891, 0.07394483817, 2.449384416, -0.06270726266, 2.379920891, 2.449384416;

  if(false)
    std::cerr << std::setprecision(10) << a.phase_function_moment_each_layer(13000).value()
      (Range(0,1), 0, Range::all()) << "\n";
  BOOST_CHECK_MATRIX_CLOSE_TOL(a.phase_function_moment_each_layer(13000).value()
                               (Range(0,1), 0, Range::all()), pf_expect, 1e-5);
}

BOOST_AUTO_TEST_SUITE_END()


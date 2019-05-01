#include "absorber_vmr_level_scaled2.h"
#include "met_data_fixture.h"
#include "unit_test_support.h"

using namespace FullPhysics;
using namespace blitz;

BOOST_FIXTURE_TEST_SUITE(absorber_vmr_level_scaled2, MetDataFixture)

BOOST_AUTO_TEST_CASE(basic)
{
  Array<double, 1> grid_expect(19);
  Array<bool, 1> scale_flag(1);
  grid_expect = 6.52839575e-06, 4.74787247e-06, 4.92904187e-06, 7.15214080e-06, 3.08558918e-05, 
    7.61965519e-05, 1.16499294e-04, 1.36102208e-04, 1.84245383e-04, 2.64503160e-04, 3.69431451e-04,
    4.66057518e-04, 5.48717085e-04, 5.73098800e-04, 6.15544205e-04, 7.41377218e-04, 
    9.09569234e-04, 1.30658765e-03, 1.46738835e-03;
  scale_flag(0) = true;

  AbsorberVmrLevelScaled2 avmr(pressure, grid_expect, 1.0, scale_flag, "H2O");

  for(int i = 0; i < grid_expect.rows(); ++i)
    BOOST_CHECK_CLOSE(avmr.volume_mixing_ratio(pressure->pressure_grid()(i).value).value(), grid_expect(i), 1e-3);
}

BOOST_AUTO_TEST_SUITE_END()

/*
 AbsorberVmrLevelScaled2::AbsorberVmrLevelScaled2(const boost::shared_ptr<Pressure>& Press,
                                                 const blitz::Array<double, 1>& Vmr_profile,
                                                 double Scale,
                                                 const blitz::Array<bool, 1>& Scale_flag,
                                                 const std::string& Gas_name)
 */

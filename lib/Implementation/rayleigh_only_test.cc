#include "l2_fp_configuration_lua.h"
#include "radiative_transfer.h"
#include "output_hdf.h"
#include "unit_test_support.h"

using namespace FullPhysics;
using namespace blitz;
BOOST_FIXTURE_TEST_SUITE(rayleigh_only, GlobalFixture)

BOOST_AUTO_TEST_CASE(radiance_and_jacobian)
{
  // Compare a Rayleigh only radiative transfer with one that has
  // aerosols in it, but all the extinction coefficients set to zero.
  L2FpConfigurationLua c_rayleigh(test_data_dir() + "lua/config_rayleigh_only.lua");
  const RadiativeTransfer& rt_rayleigh = *c_rayleigh.lua_state().globals()["config"]["rt"].value_ptr<RadiativeTransfer>();

  Array<double, 1> wn_arr(10);
  wn_arr = 12930, 12940, 12950, 12960, 12970, 12980, 12990, 13000, 13010, 13020;

  ArrayAd<double, 1> refl_rayleigh = rt_rayleigh.reflectance(wn_arr, 0).spectral_range().data_ad();

  IfstreamCs expt_data(test_data_dir() + "expected/rayleigh_only/reflectance");
  Array<double, 1> refl_val_expt;
  Array<double, 2> refl_jac_expt;
  expt_data >> refl_val_expt >> refl_jac_expt;

  BOOST_CHECK_MATRIX_CLOSE(refl_val_expt, refl_rayleigh.value());

  Array<double, 2> jac_rayleigh = refl_rayleigh.jacobian();

  BOOST_CHECK_MATRIX_CLOSE(refl_jac_expt, jac_rayleigh);
}

BOOST_AUTO_TEST_SUITE_END()

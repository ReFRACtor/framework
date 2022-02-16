#include "serialized_configuration_fixture.h"
#include "radiative_transfer.h"
#include "output_hdf.h"
#include "unit_test_support.h"

using namespace FullPhysics;
using namespace blitz;

class RayleighOnlyConfigurationFixture : public SerializedConfigurationFixture {
public:
    RayleighOnlyConfigurationFixture() : SerializedConfigurationFixture("rayleigh_only_example_config.bin.gz") { ; }
};

BOOST_FIXTURE_TEST_SUITE(rayleigh_only, RayleighOnlyConfigurationFixture)

BOOST_AUTO_TEST_CASE(radiance_and_jacobian)
{
  Array<double, 1> wn_arr(10);
  wn_arr = 12930, 12940, 12950, 12960, 12970, 12980, 12990, 13000, 13010, 13020;

  ArrayAd<double, 1> refl_rayleigh = config_rt->reflectance(wn_arr, 0).spectral_range().data_ad();

  if (false) {
    // Regenerate expected values
    std::ofstream out(test_data_dir() + "expected/rayleigh_only/reflectance");
    out << std::setprecision(20) << std::scientific;
    out << "# reflectance value" << std::endl << refl_rayleigh.value() << std::endl;
    out << "# reflectance jacobian" << std::endl << refl_rayleigh.jacobian() << std::endl;
  }

  IfstreamCs expt_data(test_data_dir() + "expected/rayleigh_only/reflectance");
  Array<double, 1> refl_val_expt;
  Array<double, 2> refl_jac_expt;
  expt_data >> refl_val_expt >> refl_jac_expt;

  BOOST_CHECK_MATRIX_CLOSE(refl_val_expt, refl_rayleigh.value());

  Array<double, 2> jac_rayleigh = refl_rayleigh.jacobian();

  BOOST_CHECK_MATRIX_CLOSE(refl_jac_expt, jac_rayleigh);
}

BOOST_AUTO_TEST_SUITE_END()

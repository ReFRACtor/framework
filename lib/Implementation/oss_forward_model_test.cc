#include "oss_forward_model.h"
#include "unit_test_support.h"
#include "oss_configuration_fixture.h"

using namespace FullPhysics;
using namespace blitz;

BOOST_FIXTURE_TEST_SUITE(oss_forward_constant, OssConfigurationFixture)

BOOST_AUTO_TEST_CASE(radiance)
{
  std::string sel_file = oss_data_dir() + "aura-tes-B1B2-unapod-loc-clear-23V-M12.4-v1.0.train.sel";
  std::string od_file = oss_data_dir() + "aura-tes-B1B2-unapod-loc-clear-23V-M12.4-v1.0.train.lut";
  std::string sol_file = oss_data_dir() + "newkur.dat";
  std::string fix_file = oss_data_dir() + "default.dat";
  std::string ch_sel_file = "NULL";
  OssForwardModel fm(config_atmosphere, config_absorber, config_absorber_calc_jacob, config_pressure, sel_file,
          od_file, sol_file, fix_file, ch_sel_file);
  fm.setup_grid();
  SpectralDomain oss_spec_domain = fm.spectral_domain(0);
  Array<double,1> oss_wavenumbers = oss_spec_domain.wavenumber();
  float expected_start_wavelength = 923.0;
  float expected_wavelength_step = 0.06;
  int expected_number_wavelength = 3951;
  BOOST_CHECK_EQUAL(oss_wavenumbers.rows(), expected_number_wavelength);
  for (int i = 0; i < expected_number_wavelength; i++) {
      BOOST_CHECK_CLOSE(oss_wavenumbers(i), expected_start_wavelength + (i * expected_wavelength_step), 1e-5);
  }
  Spectrum radiance = fm.radiance(0);
  // TODO: Add checks
}

BOOST_AUTO_TEST_SUITE_END()

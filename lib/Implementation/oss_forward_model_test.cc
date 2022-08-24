#include "oss_forward_model.h"
#include "unit_test_support.h"
#include "oss_configuration_fixture.h"

using namespace FullPhysics;
using namespace blitz;

BOOST_FIXTURE_TEST_SUITE(oss_forward_model, OssConfigurationFixture, * boost::unit_test::disabled())

BOOST_AUTO_TEST_CASE(radiance)
{
  std::string sel_file = oss_data_dir() + "aura-tes-B1B2-unapod-loc-clear-23V-M12.4-v1.0.train.sel";
  std::string od_file = oss_data_dir() + "aura-tes-B1B2-unapod-loc-clear-23V-M12.4-v1.0.train.lut";
  std::string sol_file = oss_data_dir() + "newkur.dat";
  std::string fix_file = oss_data_dir() + "default.dat";
  std::string ch_sel_file = "NULL";
  OssForwardModel fm(config_vmr, config_pressure, config_temperature, config_skin_temperature,
          config_ground, config_obs_zen_ang, config_sol_zen_ang, config_lat, config_surf_alt,
          config_lambertian, sel_file, od_file, sol_file, fix_file, ch_sel_file);

  boost::shared_ptr<OssRetrievalFlags> retrieval_flags = boost::make_shared<OssRetrievalFlags>(retrieval_temperature_levels,
          retrieval_skin_temperature_flag, retrieval_gas_levels, retrieval_emissivity_flags,
          retrieval_reflectivity_flags);
  fm.setup_retrieval(retrieval_flags);

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
  IfstreamCs expected(test_data_dir() + "expected/oss_interface/radiance_and_jacobian");

  Array<double, 1> rad_expect;
  expected >> rad_expect;
  BOOST_CHECK_MATRIX_CLOSE_TOL(radiance.spectral_range().data(), rad_expect, 1e-5);

  Array<double, 2> xk_temp_expect;
  expected >> xk_temp_expect;
  BOOST_CHECK_MATRIX_CLOSE_TOL(fm.cached_outputs->xk_temp.value, xk_temp_expect, 1e-5);

  Array<double, 1> xk_tskin_expect;
  expected >> xk_tskin_expect;
  BOOST_CHECK_MATRIX_CLOSE_TOL(fm.cached_outputs->xk_tskin.value, xk_tskin_expect, 1e-5);

  Array<double, 3> xk_out_gas_expect;
  expected >> xk_out_gas_expect;
  BOOST_CHECK_MATRIX_CLOSE_TOL(fm.cached_outputs->xk_out_gas.value, xk_out_gas_expect, 1e-5);

  Array<double, 2> xk_em_expect;
  expected >> xk_em_expect;
  BOOST_CHECK_MATRIX_CLOSE_TOL(fm.cached_outputs->xk_em.value, xk_em_expect, 1e-5);

  Array<double, 2> xk_rf_expect;
  expected >> xk_rf_expect;
  BOOST_CHECK_MATRIX_CLOSE_TOL(fm.cached_outputs->xk_rf.value, xk_rf_expect, 1e-5);

  Array<double, 1> xk_cldln_pres_expect;
  expected >> xk_cldln_pres_expect;
  BOOST_CHECK_MATRIX_CLOSE_TOL(fm.cached_outputs->xk_cldln_pres.value, xk_cldln_pres_expect, 1e-5);

  Array<double, 2> xk_cldln_ext_expect;
  expected >> xk_cldln_ext_expect;
  BOOST_CHECK_MATRIX_CLOSE_TOL(fm.cached_outputs->xk_cldln_ext.value, xk_cldln_ext_expect, 1e-5);

}

BOOST_AUTO_TEST_SUITE_END()

#include "spectrum_sampling_fixed_spacing.h"
#include "simple_fixed_spectrum_sampling.h"
#include "unit_test_support.h"
#include "configuration_fixture.h"

using namespace FullPhysics;
using namespace blitz;

BOOST_FIXTURE_TEST_SUITE(spectrum_sampling_fixed_spacing, ConfigurationFixture)

BOOST_AUTO_TEST_CASE(basic)
{
  blitz::Array<double, 1> spec_spac_val(3);
  spec_spac_val = 0.01;
  ArrayWithUnit<double, 1> spec_spac_awu(spec_spac_val, units::inv_cm);
  SpectrumSamplingFixedSpacing ssamp(spec_spac_awu);

  // Hardcoded results that we expect;
  SimpleFixedSpectrumSampling sexpect(12955.55, 13179.13, 0.01,
                                  6177.14, 6262.15, 0.01,
                                  4805.02, 4886.97, 0.01);

  // For debugging //config_spectral_window, 
  if(false) {
    for(int i = 0; i < 3; i++) {
      Array<double,1> wn(ssamp.spectral_domain(i, lowres_grid(i), high_res_extension(i)).wavenumber());
      std::cerr << std::setprecision(8) 
        << wn(wn.rows()-1) << ", " << wn(0) << std::endl;
    }
  }
  for(int i = 0; i < 3; ++i) {
    Array<double,1> expect_wl;
    expect_wl.reference(sexpect.spectral_domain(i, lowres_grid(i), high_res_extension(i)).wavelength());
    expect_wl.reverseSelf(firstDim);

    BOOST_CHECK_MATRIX_CLOSE_TOL
      (ssamp.spectral_domain(i, lowres_grid(i), high_res_extension(i)).wavelength(),
       expect_wl, 1e-6);
  }
}

BOOST_AUTO_TEST_SUITE_END()

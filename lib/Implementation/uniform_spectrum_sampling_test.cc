#include "uniform_spectrum_sampling.h"
#include "unit_test_support.h"

using namespace FullPhysics;
using namespace blitz;

BOOST_FIXTURE_TEST_SUITE(uniform_spectrum_sampling, GlobalFixture)

BOOST_AUTO_TEST_CASE(basic)
{
    Array<double, 1> spacing(1);
    spacing = 0.5;

    UniformSpectrumSampling uniform_samp(ArrayWithUnit<double, 1>(spacing, Unit("cm^-1")));

    Array<double, 1> low_res(5);
    low_res = 1, 2, 3, 4, 5;
    SpectralDomain low_res_sd(low_res, Unit("cm^-1"));

    DoubleWithUnit ils_half_width(0.5, Unit("cm^-1"));

    Array<double, 1> expect(11);
    expect = 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5, 5.5;
 
    BOOST_CHECK_EQUAL(uniform_samp.number_spectrometer(), 1);
    BOOST_CHECK_MATRIX_CLOSE(uniform_samp.spectral_domain(0, low_res_sd, ils_half_width).wavenumber(), expect);
}

BOOST_AUTO_TEST_SUITE_END()

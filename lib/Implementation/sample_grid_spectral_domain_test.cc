#include "sample_grid_spectral_domain.h"
#include "unit_test_support.h"
using namespace FullPhysics;
using namespace blitz;

BOOST_FIXTURE_TEST_SUITE(sample_grid_spectral_domain, GlobalFixture)

BOOST_AUTO_TEST_CASE(basic)
{
    /* TODO: Better test */
    Array<double, 1> sample_wavenumbers(4);
    // 200nm, 400m, 600nm, 800nm
    sample_wavenumbers = 50000.0, 50000.0 / 2.0, 50000.0 / 3.0, 50000.0 / 4.0;
    SpectralDomain spec_domain(sample_wavenumbers, units::inv_cm);
    SampleGridSpectralDomain samp_grid(spec_domain, "Test band", 4, false);
}

BOOST_AUTO_TEST_SUITE_END()

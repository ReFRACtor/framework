#include "spectral_domain.h"
#include "fp_exception.h"
#include "unit_test_support.h"

using namespace FullPhysics;
using namespace blitz;
BOOST_FIXTURE_TEST_SUITE(spectral_domain, GlobalFixture)

BOOST_AUTO_TEST_CASE(basic_test)
{
  Array<double, 1> wavenumber(5);
  wavenumber = 1000, 1000.5, 1001, 1001.5, 1002;
  Array<double, 1> wavelength(5);
  wavelength = 1e4 / wavenumber;
  SpectralDomain wn(wavenumber);
  SpectralDomain wl(wavelength, units::micron);
  BOOST_CHECK_MATRIX_CLOSE(wn.data(), wavenumber);
  BOOST_CHECK_EQUAL(wn.units().name(), "cm^-1");
  BOOST_CHECK_EQUAL(wn.type_preference(), SpectralDomain::PREFER_WAVENUMBER);
  BOOST_CHECK_MATRIX_CLOSE(wn.wavenumber(), wavenumber);
  BOOST_CHECK_MATRIX_CLOSE(wn.wavelength(), wavelength);
  BOOST_CHECK_MATRIX_CLOSE(wn.wavenumber(1 / units::m), wavenumber * 100.0);
  BOOST_CHECK_MATRIX_CLOSE(wn.wavelength(units::nm), wavelength * 1e3);
  BOOST_CHECK_MATRIX_CLOSE(wl.data(), wavelength);
  BOOST_CHECK_EQUAL(wl.units().name(), "micron");
  BOOST_CHECK_EQUAL(wl.type_preference(), SpectralDomain::PREFER_WAVELENGTH);
  BOOST_CHECK_MATRIX_CLOSE(wl.wavenumber(), wavenumber);
  BOOST_CHECK_MATRIX_CLOSE(wl.wavelength(), wavelength);
  BOOST_CHECK_MATRIX_CLOSE(wl.wavenumber(1 / units::m), wavenumber * 100.0);
  BOOST_CHECK_MATRIX_CLOSE(wl.wavelength(units::nm), wavelength * 1e3);
}

BOOST_AUTO_TEST_CASE(padding)
{
    Array<double, 1> expt_grid_fwd(90);
    double nm = 300;
    double spacing = 0.1;
    int num_jac = 5;

    for(int idx = 0; idx < expt_grid_fwd.rows(); idx++) {
        expt_grid_fwd(idx) = nm;
        nm += spacing;
    }

    Array<double, 2> expt_jac(90, num_jac);
    expt_jac(Range::all(), Range::all()) = 0.0;
    expt_jac(Range(30, 59), Range::all()) = 1.0;

    // Check increasing order
    ArrayAd<double, 1> inp_grid_fwd( expt_grid_fwd(Range(30, 59)), num_jac );
    inp_grid_fwd.jacobian() = 1.0;

    SpectralDomain inp_sd_fwd = SpectralDomain(inp_grid_fwd, units::nm);

    SpectralDomain padded_fwd = inp_sd_fwd.add_padding(DoubleWithUnit(expt_grid_fwd(30) - expt_grid_fwd(0), units::nm));

    BOOST_CHECK_MATRIX_CLOSE(expt_grid_fwd, padded_fwd.data());
    BOOST_CHECK_MATRIX_CLOSE(expt_jac, padded_fwd.data_ad().jacobian());

    // Check decreasing order
    Array<double, 1> expt_grid_rev = expt_grid_fwd.reverse(firstDim);

    ArrayAd<double, 1> inp_grid_rev( expt_grid_rev(Range(30, 59)), num_jac );
    inp_grid_rev.jacobian() = 1.0;

    SpectralDomain inp_sd_rev = SpectralDomain(inp_grid_rev, units::nm);

    SpectralDomain padded_sd = inp_sd_rev.add_padding(DoubleWithUnit(expt_grid_rev(0) - expt_grid_rev(30), units::nm));

    BOOST_CHECK_MATRIX_CLOSE(expt_grid_rev, padded_sd.data());
    BOOST_CHECK_MATRIX_CLOSE(expt_jac, padded_fwd.data_ad().jacobian());

}

BOOST_AUTO_TEST_SUITE_END()

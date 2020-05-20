#include "dispersion_polynomial.h"
#include "spectral_window_range.h"
#include <boost/foreach.hpp>
#include "unit_test_support.h"
#include "configuration_fixture.h"

using namespace FullPhysics;
using namespace blitz;

BOOST_FIXTURE_TEST_SUITE(spectral_window_range, ConfigurationFixture)

BOOST_AUTO_TEST_CASE(grid_indexes)
{
    std::vector<int> pix =
        config_spectral_window->grid_indexes(config_instrument->pixel_spectral_domain(0), 0);
    BOOST_CHECK_EQUAL((int) pix.size(), 836);
    BOOST_CHECK_EQUAL(pix.front(), 89);
    BOOST_CHECK_EQUAL(pix.back(), 924);

    pix =
        config_spectral_window->grid_indexes(config_instrument->pixel_spectral_domain(1), 1);
    BOOST_CHECK_EQUAL((int) pix.size(), 653);
    BOOST_CHECK_EQUAL(pix.front(), 210);
    BOOST_CHECK_EQUAL(pix.back(), 862);

    pix =
        config_spectral_window->grid_indexes(config_instrument->pixel_spectral_domain(2), 2);
    BOOST_CHECK_EQUAL((int) pix.size(), 812);
    BOOST_CHECK_EQUAL(pix.front(), 99);
    BOOST_CHECK_EQUAL(pix.back(), 910);
}

BOOST_AUTO_TEST_CASE(apply)
{
    SpectralDomain sdall = config_instrument->pixel_spectral_domain(0);
    SpectralDomain sd = config_spectral_window->apply(sdall, 0);
    BOOST_CHECK_CLOSE(FullPhysics::conversion(sdall.units(), sd.units()), 1.0, 1e-8);
    BOOST_CHECK_MATRIX_CLOSE(sd.data(), sdall.data()(Range(89, 924)));
}

BOOST_AUTO_TEST_CASE(apply_multi)
{
    HdfFile h(test_data_dir() + "/in/spectral_window_range/l2_multimicrowindow.h5");
    SpectralWindowRange swin(h.read_field_with_unit<double, 3>
                             ("Spectral_Window/microwindow"));
    std::vector<int> pix =
        swin.grid_indexes(config_instrument->pixel_spectral_domain(0), 0);
    BOOST_CHECK_EQUAL((int) pix.size(), 880);
    BOOST_CHECK_EQUAL(pix.front(), 29);
    BOOST_CHECK_EQUAL(pix.back(), 984);
    pix = swin.grid_indexes(config_instrument->pixel_spectral_domain(1), 1);
    BOOST_CHECK_EQUAL((int) pix.size(), 835);
    BOOST_CHECK_EQUAL(pix.front(), 7);
    BOOST_CHECK_EQUAL(pix.back(), 1015);
    pix = swin.grid_indexes(config_instrument->pixel_spectral_domain(2), 2);
    BOOST_CHECK_EQUAL((int) pix.size(), 732);
    BOOST_CHECK_EQUAL(pix.front(), 0);
    BOOST_CHECK_EQUAL(pix.back(), 887);
}

BOOST_AUTO_TEST_CASE(wavelength_file)
{
    HdfFile h(test_data_dir() + "/in/spectral_window_range/l2_ocowin.h5");
    SpectralWindowRange swin(h.read_field_with_unit<double, 3>
                             ("Spectral_Window/microwindow"));
    Array<double, 2> expect_range(3, 2);
    expect_range =
        0.755, 0.785,
        1.58, 1.65,
        2.03, 2.09;
    SpectralBound sb = swin.spectral_bound();
    for (int i = 0; i < expect_range.rows(); i++) {
        BOOST_CHECK_EQUAL(sb.lower_bound(i).value, expect_range(i, 0));
        BOOST_CHECK_EQUAL(sb.upper_bound(i).value, expect_range(i, 1));
        BOOST_CHECK_EQUAL(sb.lower_bound(i).convert_wave(units::micron).value,
                          expect_range(i, 0));
        BOOST_CHECK_EQUAL(sb.upper_bound(i).convert_wave(units::micron).value,
                          expect_range(i, 1));
    }

}

BOOST_AUTO_TEST_CASE(hdf_read)
{
  Array<double, 2> expect_range(3, 2);
  expect_range =
    0.75918488928043171, 0.77146052770058471,
    1.5979769247163076, 1.6177916201995313,
    2.0476615743155038, 2.0797545539697246;

  SpectralBound sb = config_spectral_window->spectral_bound();
  for (int i = 0; i < expect_range.rows(); i++) {
    BOOST_REQUIRE(sb.lower_bound(i).units.name() == "micron");
    BOOST_CHECK_CLOSE(sb.lower_bound(i).value, expect_range(i, 0), 1e-6);
    BOOST_CHECK_CLOSE(sb.upper_bound(i).value, expect_range(i, 1), 1e-6);
  }
}

BOOST_AUTO_TEST_CASE(bounds_ordering)
{
  SpectralBound sb = config_spectral_window->spectral_bound();
  for(int i = 0; i < 3; i++) {
    BOOST_REQUIRE(sb.lower_bound(i) < sb.upper_bound(i));
    BOOST_REQUIRE(sb.lower_bound(i, units::inv_cm) < 
                  sb.upper_bound(i, units::inv_cm));
    }
}

BOOST_AUTO_TEST_CASE(bad_sample_mask_constructor_1)
{
    // Use the typical OCO all band ranges
    Array<double, 3> win_ranges(3, 1, 2);
    win_ranges =
        0.755, 0.785,
        1.58, 1.65,
        2.03, 2.09;

    // Use a double to replicate what would come from L1B file in snr_coeff,
    // No reason we couldn't just use a bool mask
    Array<bool, 2> bad_sample_mask(3, 1016);
    bad_sample_mask = false;

    // Mark some bad samples
    bad_sample_mask(0, Range(0,100)) = true;
    bad_sample_mask(1, Range(10,20)) = true;
    bad_sample_mask(2, Range(0, 9)) = true;
    bad_sample_mask(2, Range(1006, 1015)) = true;

    SpectralWindowRange spec_win_range(ArrayWithUnit<double, 3>(win_ranges, Unit("um")), bad_sample_mask);

    std::vector<int> pix = spec_win_range.grid_indexes(config_instrument->pixel_spectral_domain(0), 0);
    BOOST_CHECK_EQUAL((int) pix.size(), 1016-101);
    for(int i = 0; i < (int) pix.size(); i++) {
        BOOST_CHECK_EQUAL(pix[i], i+101);
    }

    pix = spec_win_range.grid_indexes(config_instrument->pixel_spectral_domain(1), 1);
    BOOST_CHECK_EQUAL((int) pix.size(), 1016-11);
    for(int i = 0; i < (int) pix.size(); i++) {
        if(i < 10) {
            BOOST_CHECK_EQUAL(pix[i], i);
        } else {
            BOOST_CHECK_EQUAL(pix[i], i+11);
        }
    }

    pix = spec_win_range.grid_indexes(config_instrument->pixel_spectral_domain(2), 2);
    BOOST_CHECK_EQUAL((int) pix.size(), 1016-20);
    for(int i = 0; i < (int) pix.size(); i++) {
        BOOST_CHECK_EQUAL(pix[i], i+10);
    }

}

BOOST_AUTO_TEST_CASE(bad_sample_mask_constructor_2)
{
    // Use the typical OCO all band ranges
    Array<double, 3> win_ranges(3, 1, 2);
    win_ranges =
        0.755, 0.785,
        1.58, 1.65,
        2.03, 2.09;

    // Use a double to replicate what would come from L1B file in snr_coeff,
    // No reason we couldn't just use a bool mask
    Array<bool, 1> bad_sample_mask_1(1016);
    bad_sample_mask_1 = false;

    Array<bool, 1> bad_sample_mask_2(1016);
    bad_sample_mask_2 = false;

    Array<bool, 1> bad_sample_mask_3(1016);
    bad_sample_mask_3 = false;

    // Mark some bad samples
    bad_sample_mask_1(Range(0,100)) = true;
    bad_sample_mask_2(Range(10,20)) = true;
    bad_sample_mask_3(Range(0, 9)) = true;
    bad_sample_mask_3(Range(1006, 1015)) = true;

    std::vector<Array<bool, 1> > bad_sample_mask;
    bad_sample_mask.push_back(bad_sample_mask_1);
    bad_sample_mask.push_back(bad_sample_mask_2);
    bad_sample_mask.push_back(bad_sample_mask_3);

    SpectralWindowRange spec_win_range(ArrayWithUnit<double, 3>(win_ranges, Unit("um")), bad_sample_mask);

    std::vector<int> pix = spec_win_range.grid_indexes(config_instrument->pixel_spectral_domain(0), 0);
    BOOST_CHECK_EQUAL((int) pix.size(), 1016-101);
    for(int i = 0; i < (int) pix.size(); i++) {
        BOOST_CHECK_EQUAL(pix[i], i+101);
    }

    pix = spec_win_range.grid_indexes(config_instrument->pixel_spectral_domain(1), 1);
    BOOST_CHECK_EQUAL((int) pix.size(), 1016-11);
    for(int i = 0; i < (int) pix.size(); i++) {
        if(i < 10) {
            BOOST_CHECK_EQUAL(pix[i], i);
        } else {
            BOOST_CHECK_EQUAL(pix[i], i+11);
        }
    }

    pix = spec_win_range.grid_indexes(config_instrument->pixel_spectral_domain(2), 2);
    BOOST_CHECK_EQUAL((int) pix.size(), 1016-20);
    for(int i = 0; i < (int) pix.size(); i++) {
        BOOST_CHECK_EQUAL(pix[i], i+10);
    }

}


BOOST_AUTO_TEST_CASE(empty_spectral_bounds)
{
  // Test that when we have empty spectral bounds that there is no error
  // just because they are zeroed out
  Array<double, 3> win_range(3, 1, 2);
  win_range =
      -1, -1,
      0, 0,
      1, 1;

  // Create dispersion vector
  std::vector<boost::shared_ptr<SampleGrid> > spec_disp;
  Array<bool, 1> disp_flag(2);
  disp_flag = true, false;
  Array<double, 1> disp_coeff(2);

  disp_coeff = 0.757691, 1.74757e-05;
  spec_disp.push_back(boost::shared_ptr<DispersionPolynomial>(new DispersionPolynomial(disp_coeff, disp_flag, units::micron, "ABO2", 1016, true)));

  disp_coeff = 1.59071, 3.62647e-05;
  spec_disp.push_back(boost::shared_ptr<DispersionPolynomial>(new DispersionPolynomial(disp_coeff, disp_flag, units::micron, "WCO2", 1016, true)));

  disp_coeff = 2.04325, 4.69383e-05;
  spec_disp.push_back(boost::shared_ptr<DispersionPolynomial>(new DispersionPolynomial(disp_coeff, disp_flag, units::micron, "SCO2", 1016, true)));


  SpectralWindowRange spec_win = SpectralWindowRange(ArrayWithUnit<double, 3>(win_range, units::sample_index));
  spec_win.dispersion(spec_disp);

  SpectralBound sb = spec_win.spectral_bound();
  for (int i = 0; i < win_range.rows(); i++) {
    BOOST_REQUIRE(sb.lower_bound(i).units.name() == "micron");
    BOOST_CHECK_EQUAL(sb.lower_bound(i).value - sb.upper_bound(i).value, 0);
  }

}

BOOST_AUTO_TEST_SUITE_END()

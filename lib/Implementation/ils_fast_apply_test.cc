#include "ils_fast_apply.h"
#include "hdf_file.h"
#include "dispersion_polynomial.h"
#include "unit_test_support.h"

using namespace FullPhysics;
using namespace blitz;

BOOST_FIXTURE_TEST_SUITE(ils_fast_apply, GlobalFixture)

BOOST_AUTO_TEST_CASE(basic)
{
    int num_pixels = 1016;

    std::vector<int> pixel_list;
    for(int i = 0; i < num_pixels; i++) {
        pixel_list.push_back(i);
    }

    Array<bool, 1> flag(6);
    flag = false;
  
    Array<double, 1> coeff(6);
    coeff =  7.57646001e-01,   1.74935138e-05,  -2.73614810e-09,
        -7.51159779e-14,   1.32166257e-16,  -7.27469045e-20;
  
    boost::shared_ptr<DispersionPolynomial>
      dispersion(new DispersionPolynomial(coeff, flag, units::micron, "A-Band", num_pixels, true));
  
    HdfFile inp_vals(test_data_dir() + "expected/ils_fast_apply/fft_svd_ils_values.h5");
  
    auto scaled_uh_isrf = inp_vals.read_field<double, 2>("scaled_uh_isrf");
    auto svh_isrf_fft_real = inp_vals.read_field<double, 2>("svh_isrf_fft_real");
    auto svh_isrf_fft_imag = inp_vals.read_field<double, 2>("svh_isrf_fft_imag");
    auto extract_indices = inp_vals.read_field<int, 1>("center_freq_indices");
  
    auto ils_half_width = DoubleWithUnit(4.09e-04, "um");
  
    auto ils = IlsFastApply(scaled_uh_isrf,
                         svh_isrf_fft_real,
                         svh_isrf_fft_imag,
                         extract_indices, 
                         dispersion, ils_half_width, "A-Band", "abo2");

    auto high_res_grid = inp_vals.read_field<double, 1>("high_res_grid");
    auto high_res_radiance = inp_vals.read_field<double, 1>("high_res_radiance");

    ArrayAd<double, 1> high_res_rad_ad(high_res_radiance.rows(), 0);
    high_res_rad_ad.value() = high_res_radiance;

    auto calc_conv_rad = ils.apply_ils(high_res_grid, high_res_rad_ad, pixel_list);

    auto expected_conv_rad = inp_vals.read_field<double, 1>("expected_convolved_radiance");

    BOOST_CHECK_MATRIX_CLOSE_TOL(expected_conv_rad, calc_conv_rad.value(), 1e-7);
}

BOOST_AUTO_TEST_SUITE_END()

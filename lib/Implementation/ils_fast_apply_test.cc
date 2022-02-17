#include "ils_fast_apply.h"
#include "unit_test_support.h"
#include "hdf_file.h"
#include "sample_grid_spectral_domain.h"

using namespace FullPhysics;
using namespace blitz;

BOOST_FIXTURE_TEST_SUITE(ils_fast_apply, GlobalFixture)

BOOST_AUTO_TEST_CASE(omi_muses_case)
{
    // Extract expected values and inputs
    HdfFile inp_vals(test_data_dir() + "expected/ils_fast_apply/omi_muses_case.h5");
  
    Array<double, 2> svh_isrf_fft_real = inp_vals.read_field<double, 2>("ILS_Parameters/svh_isrf_fft_real");
    Array<double, 2> svh_isrf_fft_imag = inp_vals.read_field<double, 2>("ILS_Parameters/svh_isrf_fft_imag");
    Array<double, 2> scaled_uh_isrf = inp_vals.read_field<double, 2>("ILS_Parameters/scaled_uh_isrf");
    Array<int, 1> extract_indices = inp_vals.read_field<int, 1>("ILS_Parameters/extract_indices");

    Array<double, 1> mono_grid = inp_vals.read_field<double, 1>("Monochromatic/grid");
    Array<double, 1> mono_rad = inp_vals.read_field<double, 1>("Monochromatic/radiance");
    Array<double, 2> mono_jac = inp_vals.read_field<double, 2>("Monochromatic/jacobian");

    Array<double, 1> conv_expt_grid = inp_vals.read_field<double, 1>("Convolved/grid");
    Array<double, 1> conv_expt_rad = inp_vals.read_field<double, 1>("Convolved/radiance");
    Array<double, 2> conv_expt_jac = inp_vals.read_field<double, 2>("Convolved/jacobian");

    // Some hard coded values
    std::string band_name = "Band_1";
    DoubleWithUnit high_res_ext = DoubleWithUnit(0.315, units::nm);

    // Set up SampleGrid object
    SpectralDomain conv_sd(conv_expt_grid, units::nm);
    boost::shared_ptr<SampleGridSpectralDomain> samp_grid(new SampleGridSpectralDomain(conv_sd, band_name));

    // Create IlsFastApply object
    IlsFastApply ils = IlsFastApply(scaled_uh_isrf, svh_isrf_fft_real, svh_isrf_fft_imag, extract_indices, samp_grid, high_res_ext, band_name, band_name);

    // Create a pixel list that references all samples on the convolved grid
    std::vector<int> pixel_list;
    for(int i = 0; i < conv_expt_grid.rows(); i++) {
        pixel_list.push_back(i);
    }

    // Convert mono grid units to wavenumbers to meet interface requirements
    Array<double, 1> mono_wn = ArrayWithUnit<double, 1>(mono_grid, units::nm).convert_wave(units::inv_cm).value;

    // Apply ILS convolution for no jacobian case
    Array<double, 1> conv_calc_rad = ils.apply_ils(mono_wn, mono_rad, pixel_list);

    BOOST_CHECK_MATRIX_CLOSE_TOL(conv_expt_rad, conv_calc_rad, 1e-10);

    // Apply ILS convolution for jacobian case
    ArrayAd<double, 1> mono_rad_jac(mono_rad, mono_jac);

    ArrayAd<double, 1> conv_calc_rad_jac = ils.apply_ils(mono_wn, mono_rad_jac, pixel_list);

    BOOST_CHECK_MATRIX_CLOSE_TOL(conv_expt_rad, conv_calc_rad_jac.value(), 1e-10);
    BOOST_CHECK_MATRIX_CLOSE_TOL(conv_expt_jac, conv_calc_rad_jac.jacobian(), 1e-10);
}

BOOST_AUTO_TEST_SUITE_END()

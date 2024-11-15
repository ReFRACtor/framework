#include "ils_fast_apply.h"
#include "ostream_pad.h"

using namespace FullPhysics;
using namespace blitz;

//-----------------------------------------------------------------------
// scaled_uh_isrf: M x N
// svh_isrf_fft P x N
// M = Number of pixels, N = singular values, P = FFT size / spectal dim
//-----------------------------------------------------------------------

IlsFastApply::IlsFastApply(const blitz::Array<double, 2>& Scaled_uh_isrf,
                           const blitz::Array<double, 2>& Svh_isrf_fft_real,
                           const blitz::Array<double, 2>& Svh_isrf_fft_imag,
                           const blitz::Array<int, 1>& Extract_indices,
                           const boost::shared_ptr<SampleGrid>& Sample_grid,
                           const DoubleWithUnit& High_res_extension,
                           const std::string& Band_name, const std::string& Hdf_band_name)
:   sample_grid(Sample_grid), high_res_extension_(High_res_extension),
    band_name_(Band_name), hdf_band_name_(Hdf_band_name)
{
    // Some sanity checking on input array sizes
    if (Scaled_uh_isrf.rows() != Extract_indices.rows()) {
        Exception err;
        err << "scaled_uh_isrf rows: " << Scaled_uh_isrf.rows() 
            << " must match extract_indices rows: " << Extract_indices.rows();
        throw err;
    }

    if (Scaled_uh_isrf.cols() != Svh_isrf_fft_real.cols()) {
        Exception err;
        err << "scaled_uh_isrf cols: " << Scaled_uh_isrf.cols() 
            << " must match svh_isrf_fft cols: " << Svh_isrf_fft_real.cols();
        throw err;
    }

    if (Svh_isrf_fft_real.rows() != Svh_isrf_fft_imag.rows()) {
        Exception err;
        err << "svh_isrf_fft_real rows: " << Svh_isrf_fft_real.rows()
            << " must match svh_isrf_fft_imag rows: " << Svh_isrf_fft_imag.rows();
        throw err;
    }

    if (Svh_isrf_fft_real.cols() != Svh_isrf_fft_imag.cols()) {
        Exception err;
        err << "svh_isrf_fft_real cols: " << Svh_isrf_fft_real.cols()
            << " must match svh_isrf_fft_imag cols: " << Svh_isrf_fft_imag.cols();
        throw err;
    }

    scaled_uh_isrf.reference(Scaled_uh_isrf);

    svh_isrf_fft.resize(Svh_isrf_fft_real.shape());
    real(svh_isrf_fft) = Svh_isrf_fft_real;
    imag(svh_isrf_fft) = Svh_isrf_fft_imag;

    extract_indices.reference(Extract_indices);

    fft_size = svh_isrf_fft.rows();
    fft_wavetable = gsl_fft_complex_wavetable_alloc(fft_size);
    fft_workspace = gsl_fft_complex_workspace_alloc(fft_size);
}

IlsFastApply::~IlsFastApply()
{
    gsl_fft_complex_wavetable_free(fft_wavetable);
    gsl_fft_complex_workspace_free(fft_workspace);
}

//-----------------------------------------------------------------------
// Apply the ILS using SVD components
//-----------------------------------------------------------------------

blitz::Array<double, 1> IlsFastApply::apply_ils
(const blitz::Array<double, 1>& UNUSED(high_resolution_wave_number),
 const blitz::Array<double, 1>& high_resolution_radiance,
 const std::vector<int>& pixel_list) const
{
    // Check that the number of extract indexes is the same length as the pixel list
    if ((int) pixel_list.size() != extract_indices.rows()) {
        Exception err;
        err << "The number of samples in the pixel_list: " << pixel_list.size() << " "
            << "must match the number of extract_indices: " << extract_indices.rows();
        throw err;
    }

    // Output array
    Array<double, 1> conv_rad(pixel_list.size());
    conv_rad = 0.0;

    // Allocate array for use in FFTs
    Array<complex<double>, 1> spectrum_fft(fft_size);

    // Compute the FFT of the high resolution radiance data
    real(spectrum_fft)(Range(0, high_resolution_radiance.rows()-1)) = high_resolution_radiance;
    real(spectrum_fft)(Range(high_resolution_radiance.rows(), fft_size-1)) = 0.0;
    imag(spectrum_fft) = 0.0;

    int res_fwd = gsl_fft_complex_forward(reinterpret_cast<double *>(spectrum_fft.dataFirst()), 1, fft_size, fft_wavetable, fft_workspace);
    if (res_fwd != GSL_SUCCESS) {
        Exception err;
        err << "Error calling gsl_fft_complex_forward: " << gsl_strerror(res_fwd);
        throw err;
    }

    // Apply SVD factors to high resolution specta
    // num_s = number of singular values
    int num_s = svh_isrf_fft.cols();
    Array<complex<double>, 1> inv_fft(fft_size);
    for(int s_idx = 0; s_idx < num_s; s_idx++) {
        inv_fft = svh_isrf_fft(Range::all(), s_idx) * spectrum_fft;

        // Do inverse FFT to get convolved radiance
        int res_inv = gsl_fft_complex_inverse(reinterpret_cast<double *>(inv_fft.dataFirst()), 1, fft_size, fft_wavetable, fft_workspace);
        if (res_inv != GSL_SUCCESS) {
            Exception err;
            err << "Error calling gsl_fft_complex_inverse: " << gsl_strerror(res_inv);
            throw err;
        }

        for(int list_idx = 0; list_idx < (int) pixel_list.size(); list_idx++) {
            int pix_idx = pixel_list[list_idx];
            int hr_idx = extract_indices(pix_idx);
            conv_rad(list_idx) += scaled_uh_isrf(pix_idx, s_idx) * inv_fft[0](hr_idx);
        }
    }

    return conv_rad;
}

//-----------------------------------------------------------------------
// Apply the ILS to the radiance and jacobians
//-----------------------------------------------------------------------
                      
ArrayAd<double, 1> IlsFastApply::apply_ils(const blitz::Array<double, 1>& high_resolution_wave_number,
                                           const ArrayAd<double, 1>& high_resolution_radiance,
                                           const std::vector<int>& pixel_list) const
{
    ArrayAd<double, 1> conv_rad(pixel_list.size(), high_resolution_radiance.number_variable());

    conv_rad.value() = apply_ils(high_resolution_wave_number, high_resolution_radiance.value(), pixel_list);

    for(int jac_idx = 0; jac_idx < high_resolution_radiance.number_variable(); jac_idx++) {
        conv_rad.jacobian()(Range::all(), jac_idx) = apply_ils(high_resolution_wave_number, high_resolution_radiance.jacobian()(Range::all(), jac_idx), pixel_list);
    }

    return conv_rad;
}

void IlsFastApply::print(std::ostream& os) const
{
    os << "IlsFastApply\n";
    OstreamPad opad(os, "  ");
    opad << "Band Name: " << band_name() << "\n"
         << "Num Singlar Values: " << svh_isrf_fft.cols() << "\n"
         << "Sample Grid: " << *sample_grid << "\n";
    opad.strict_sync();
}

boost::shared_ptr<Ils> IlsFastApply::clone() const
{
    return boost::shared_ptr<Ils>(new IlsFastApply(scaled_uh_isrf, real(svh_isrf_fft), imag(svh_isrf_fft), 
                                                   extract_indices, sample_grid, high_res_extension_, band_name_, hdf_band_name_)); 
}

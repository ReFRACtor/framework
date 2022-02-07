#include "ils_svd_fft.h"

using namespace FullPhysics;
using namespace blitz;

//-----------------------------------------------------------------------
/// left_matrix_truncated: M x N
/// right_matrix_fourier_transforms P x N
/// M = Number of pixels, N = singular values, P = FFT size / spectal dim
//-----------------------------------------------------------------------

IlsSvdFft::IlsSvdFft(const blitz::Array<double, 2>& Left_matrix_truncated,
                     const blitz::Array<double, 2>& Right_matrix_fourier_transforms_real,
                     const blitz::Array<double, 2>& Right_matrix_fourier_transforms_imag,
                     const blitz::Array<int, 1>& Center_freq_indices,
                     const boost::shared_ptr<SampleGrid>& Sample_grid, 
                     const DoubleWithUnit& High_res_extension,
                     const std::string& Band_name, const std::string& Hdf_band_name) 
:   sample_grid(Sample_grid), high_res_extension_(High_res_extension),
    band_name_(Band_name), hdf_band_name_(Hdf_band_name)
{
    // Some sanity checking on input array sizes
    if (Left_matrix_truncated.rows() != Center_freq_indices.rows()) {
        Exception err;
        err << "left_matrix_truncated rows: " << Left_matrix_truncated.rows() 
            << " must match center_freq_indices rows: " << Center_freq_indices.rows();
        throw err;
    }

    if (Left_matrix_truncated.cols() != Right_matrix_fourier_transforms_real.cols()) {
        Exception err;
        err << "left_matrix_truncated cols: " << Left_matrix_truncated.cols() 
            << " must match right_matrix_fourier_transforms cols: " << Right_matrix_fourier_transforms_real.cols();
        throw err;
    }

    if (Right_matrix_fourier_transforms_real.rows() != Right_matrix_fourier_transforms_imag.rows()) {
        Exception err;
        err << "right_matrix_fourier_transforms_real rows: " << Right_matrix_fourier_transforms_real.rows()
            << " must match right_matrix_fourier_transforms_imag rows: " << Right_matrix_fourier_transforms_imag.rows();
        throw err;
    }

    if (Right_matrix_fourier_transforms_real.cols() != Right_matrix_fourier_transforms_imag.cols()) {
        Exception err;
        err << "right_matrix_fourier_transforms_real cols: " << Right_matrix_fourier_transforms_real.cols()
            << " must match right_matrix_fourier_transforms_imag cols: " << Right_matrix_fourier_transforms_imag.cols();
        throw err;
    }

    left_matrix_truncated.reference(Left_matrix_truncated);

    right_matrix_fourier_transforms.resize(Right_matrix_fourier_transforms_real.shape());
    real(right_matrix_fourier_transforms) = Right_matrix_fourier_transforms_real;
    imag(right_matrix_fourier_transforms) = Right_matrix_fourier_transforms_imag;

    center_freq_indices.reference(Center_freq_indices);

    fft_size = right_matrix_fourier_transforms.rows();
    fft_wavetable = gsl_fft_complex_wavetable_alloc(fft_size);
    fft_workspace = gsl_fft_complex_workspace_alloc(fft_size);
}

IlsSvdFft::~IlsSvdFft()
{
    gsl_fft_complex_wavetable_free(fft_wavetable);
    gsl_fft_complex_workspace_free(fft_workspace);
}

//-----------------------------------------------------------------------
/// 
//-----------------------------------------------------------------------

blitz::Array<double, 1> IlsSvdFft::apply_ils(const blitz::Array<double, 1>& high_resolution_wave_number,
                                             const blitz::Array<double, 1>& high_resolution_radiance,
                                             const std::vector<int>& pixel_list) const
{
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
    int num_s = right_matrix_fourier_transforms.cols();
    Array<complex<double>, 1> inv_fft(fft_size);
    for(int s_idx = 0; s_idx < num_s; s_idx++) {
        inv_fft = right_matrix_fourier_transforms(Range::all(), s_idx) * spectrum_fft;

        // Do inverse FFT to get convolved radiance
        int res_inv = gsl_fft_complex_inverse(reinterpret_cast<double *>(inv_fft.dataFirst()), 1, fft_size, fft_wavetable, fft_workspace);
        if (res_inv != GSL_SUCCESS) {
            Exception err;
            err << "Error calling gsl_fft_complex_inverse: " << gsl_strerror(res_inv);
            throw err;
        }
 
        for(int list_idx = 0; list_idx < pixel_list.size(); list_idx++) {
            int pix_idx = pixel_list[list_idx];
            int hr_idx = center_freq_indices(pix_idx);
            conv_rad(list_idx) += left_matrix_truncated(pix_idx, s_idx) * inv_fft[0](hr_idx);
        }
    }

    return conv_rad;
}

//-----------------------------------------------------------------------
/// 
//-----------------------------------------------------------------------
                      
ArrayAd<double, 1> IlsSvdFft::apply_ils(const blitz::Array<double, 1>& high_resolution_wave_number,
                                        const ArrayAd<double, 1>& high_resolution_radiance,
                                        const std::vector<int>& pixel_list) const
{
    // Output array
    ArrayAd<double, 1> conv_rad(pixel_list.size(), high_resolution_radiance.number_variable());
    conv_rad = 0.0;

    // Allocate array for use in FFTs
    ArrayAd<complex<double>, 1> spectrum_fft(fft_size, high_resolution_radiance.number_variable());

    // Compute the FFT of the high resolution radiance data
    real(spectrum_fft.value())(Range(0, high_resolution_radiance.rows()-1)) = high_resolution_radiance.value();
    real(spectrum_fft.value())(Range(high_resolution_radiance.rows(), fft_size-1)) = 0.0;
    imag(spectrum_fft.value()) = 0.0;

    if (!high_resolution_radiance.is_constant()) {
        real(spectrum_fft.jacobian())(Range(0, high_resolution_radiance.rows()-1)) = high_resolution_radiance.jacobian();
        real(spectrum_fft.jacobian())(Range(high_resolution_radiance.rows(), fft_size-1)) = 0.0;
        imag(spectrum_fft.jacobian()) = 0.0;
    }

    int res_fwd = gsl_fft_complex_forward(reinterpret_cast<double *>(spectrum_fft.value().dataFirst()), 1, fft_size, fft_wavetable, fft_workspace);
    if (res_fwd != GSL_SUCCESS) {
        Exception err;
        err << "Error calling gsl_fft_complex_forward: " << gsl_strerror(res_fwd);
        throw err;
    }
 
    // Apply SVD factors to high resolution specta
    // num_s = number of singular values
    int num_s = right_matrix_fourier_transforms.cols();
    ArrayAd<complex<double>, 1> inv_fft(fft_size, spectrum_fft.number_variable());
    for(int s_idx = 0; s_idx < num_s; s_idx++) {
        inv_fft.value() = right_matrix_fourier_transforms(Range::all(), s_idx) * spectrum_fft.value();

        // Do inverse FFT to get convolved radiance
        int res_inv = gsl_fft_complex_inverse(reinterpret_cast<double *>(inv_fft.value().dataFirst()), 1, fft_size, fft_wavetable, fft_workspace);
        if (res_inv != GSL_SUCCESS) {
            Exception err;
            err << "Error calling gsl_fft_complex_inverse: " << gsl_strerror(res_inv);
            throw err;
        }
 
        for(int list_idx = 0; list_idx < pixel_list.size(); list_idx++) {
            int pix_idx = pixel_list[list_idx];
            int hr_idx = center_freq_indices(pix_idx);
            conv_rad.value()(list_idx) += left_matrix_truncated(pix_idx, s_idx) * real(inv_fft.value())(hr_idx);
        }
    }

    return conv_rad;
}

void IlsSvdFft::print(std::ostream& os) const
{
    os << "IlsSvdFft\n";
}

boost::shared_ptr<Ils> IlsSvdFft::clone() const
{
    return boost::shared_ptr<Ils>(new IlsSvdFft(left_matrix_truncated, real(right_matrix_fourier_transforms), imag(right_matrix_fourier_transforms), 
                                                center_freq_indices, sample_grid, high_res_extension_, band_name_, hdf_band_name_)); 
}

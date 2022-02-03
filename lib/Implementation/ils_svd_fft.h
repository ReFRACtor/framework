#ifndef ILS_SVD_FFT_H
#define ILS_SVD_FFT_H
#include "ils.h"
#include "dispersion.h"
#include <gsl/gsl_fft_complex.h>

namespace FullPhysics {
/****************************************************************//**
  This is a ILS where pre computed SVD values are used together
  with FFT calls to apply the convolution.
*******************************************************************/

class IlsSvdFft : public Ils {
public:
    //-----------------------------------------------------------------------
    /// Constructor with SVD array computed offline supplied.
    //-----------------------------------------------------------------------
    IlsSvdFft(const blitz::Array<double, 2>& Left_matrix_truncated,
              const blitz::Array<double, 2>& Right_matrix_fourier_transforms_real,
              const blitz::Array<double, 2>& Right_matrix_fourier_transforms_imag,
              const blitz::Array<int, 1>& Center_freq_indices,
              const boost::shared_ptr<Dispersion>& Disp, 
              const DoubleWithUnit& Ils_half_width,
              const std::string& Band_name, const std::string& Hdf_band_name);
  
    virtual ~IlsSvdFft();
  
    virtual blitz::Array<double, 1> apply_ils(const blitz::Array<double, 1>& high_resolution_wave_number,
                                              const blitz::Array<double, 1>& high_resolution_radiance,
                                              const std::vector<int>& pixel_list) const;
                          
    virtual ArrayAd<double, 1> apply_ils(const blitz::Array<double, 1>& high_resolution_wave_number,
                                         const ArrayAd<double, 1>& high_resolution_radiance,
                                         const std::vector<int>& pixel_list) const;
  
    virtual void print(std::ostream& os) const;
  
    virtual std::string band_name() const 
    { return band_name_; }
  
    virtual std::string hdf_band_name() const 
    { return hdf_band_name_; }
  
    virtual SpectralDomain pixel_grid() const
    { return disp->pixel_grid(); }
  
    virtual DoubleWithUnit ils_half_width() const
    {return ils_half_width_;}
  
    virtual void ils_half_width(const DoubleWithUnit& half_width) 
    { ils_half_width_ = half_width; }
  
    virtual boost::shared_ptr<Ils> clone() const;
  
private:
    blitz::Array<double, 2> left_matrix_truncated;
    blitz::Array<std::complex<double>, 2> right_matrix_fourier_transforms;
    blitz::Array<int, 1> center_freq_indices;
    boost::shared_ptr<Dispersion> disp;
    DoubleWithUnit ils_half_width_;
    std::string band_name_, hdf_band_name_;
    
    int fft_size;
    gsl_fft_complex_wavetable *fft_wavetable;
    gsl_fft_complex_workspace *fft_workspace;
};
}
#endif

#ifndef ILS_SVD_FFT_H
#define ILS_SVD_FFT_H
#include "ils.h"
#include "sample_grid.h"
#include <gsl/gsl_fft_complex.h>

namespace FullPhysics {
/****************************************************************//**
  This is a ILS where pre computed SVD values are used together
  with FFT calls to apply the convolution.
*******************************************************************/

class IlsFastApply : public Ils {
public:
    //-----------------------------------------------------------------------
    /// Constructor with SVD array computed offline supplied.
    //-----------------------------------------------------------------------
    IlsFastApply(const blitz::Array<double, 2>& Scaled_uh_isrf,
                 const blitz::Array<double, 2>& Svh_isrf_fft_real,
                 const blitz::Array<double, 2>& Svh_isrf_fft_imag,
                 const blitz::Array<int, 1>& Extract_indices,
                 const boost::shared_ptr<SampleGrid>& Sample_grid,
                 const DoubleWithUnit& High_res_extension,
                 const std::string& Band_name, const std::string& Hdf_band_name);
  
    virtual ~IlsFastApply();
  
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
    { return sample_grid->pixel_grid(); }
  
    virtual DoubleWithUnit high_res_extension() const
    {return high_res_extension_;}
  
    virtual void high_res_extension(const DoubleWithUnit& half_width) 
    { high_res_extension_ = half_width; }
  
    virtual boost::shared_ptr<Ils> clone() const;
  
private:
    blitz::Array<double, 2> scaled_uh_isrf;
    blitz::Array<std::complex<double>, 2> svh_isrf_fft;
    blitz::Array<int, 1> extract_indices;
    boost::shared_ptr<SampleGrid> sample_grid;
    DoubleWithUnit high_res_extension_;
    std::string band_name_, hdf_band_name_;
    
    int fft_size;
    gsl_fft_complex_wavetable *fft_wavetable;
    gsl_fft_complex_workspace *fft_workspace;
};
}
#endif

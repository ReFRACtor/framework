#ifndef ILS_SVD_FFT_H
#define ILS_SVD_FFT_H
#include "ils.h"
#include "sample_grid.h"
#include <gsl/gsl_fft_complex.h>

namespace FullPhysics {
/****************************************************************//**
  One approach to represent the Instrument Line Sphape for an
  instrument is to measure or model the response function 
  corresponding to each point of an instrument measurement grid
  or some other user-defined wavelength grid (possibly created i
  by aggregating the response functions of multiple detector pixels). 
  The result is a set of ILS functions describing the contribution
  of each wavelength's spectral radiance to that measured at the 
  center wavelength of the ILS. Each of these ILS functions
  must be multiplied with the modeled radiance spectrum. An integral 
  of this product will yield the measured radiance at the associated 
  center wavelength of the ILS. 
  
  Since the true spectral radiance is typically modeled at a high 
  resolution, this can produce a very large number of multiplications 
  to perform, which can limit the bandwidth of our ability to perform 
  retrievals. 
  
  The approach used by this module takes advantage of preprocessing that
  seeks a common structure to the ILS curves via a singular value 
  decomposition, without assumming any a priori knowledge of the ILS 
  functional form.

  This module performs a fast convolutions of the spectra with the 
  dominating singular vectors using a FFT operation. In practice, 
  this method has an improved running time by at least an order of 
  magnitude.

  The input values to this class are determined in an offline process
  that interpolates response functions to the high resolution 
  monochrmatic grid. Taking advantage of the similarity of response function
  forms a certain number of decomposition components are created
  that can them be multiplied by a FFT of the monochromatic spectra
  then scaled after an inverse FFT.
*******************************************************************/

class IlsFastApply : public Ils {
public:
    //-----------------------------------------------------------------------
    /// Construct the class using SVD values computed offline. In the variable
    /// names ISRF means instrument shape response function.
    ///
    /// scaled_uh_isrf: M x N
    ///   During preprocessing, each ILS curve is normalized. These values
    ///   capture the magnitude of the scaling so it can be undone post 
    ///   convolution.
    ///
    /// svh_isrf_fft_real, svh_isrf_fft_imag: P x N
    ///   The real and imaginary components of the FFT of the S * VH
    ///   components from an SVD decomposition algorithm.
    ///
    /// extract_indices M: 
    ///    The subset of final grid that correspond to the ILS curves, 
    ///    indexes of the applicable values to extract from the reverse FFT.
    ///
    /// Dimensions:
    /// M = Number of pixels, N = singular values, P = FFT size / spectal dim
    //-----------------------------------------------------------------------
    IlsFastApply(const blitz::Array<double, 2>& Scaled_uh_isrf,
                 const blitz::Array<double, 2>& Svh_isrf_fft_real,
                 const blitz::Array<double, 2>& Svh_isrf_fft_imag,
                 const blitz::Array<int, 1>& Extract_indices,
                 const boost::shared_ptr<SampleGrid>& Sample_grid,
                 const DoubleWithUnit& High_res_extension,
                 const std::string& Band_name, const std::string& Hdf_band_name);
  
    virtual ~IlsFastApply();

    //-----------------------------------------------------------------------
    /// Apply the ILS using offline computed FFT transformed singular
    /// value decomposition factors using this process:
    /// 1. Perform an FFT on the monochromatic radiance passed in
    /// 2. Convolve the radiance with the S*VH values in FFT space through
    ///    a simple multiplication.
    /// 3. Perform an inverse FFT on the previous product
    /// 4. Unapply ILS curve normalization 
    /// 5. Extract only those values from the inverse FFT that correspond 
    ///    to the destination grid.
    //-----------------------------------------------------------------------
    virtual blitz::Array<double, 1> apply_ils(const blitz::Array<double, 1>& high_resolution_wave_number,
                                              const blitz::Array<double, 1>& high_resolution_radiance,
                                              const std::vector<int>& pixel_list) const;

    //-----------------------------------------------------------------------
    /// Apply the ILS using offline computed FFT transformed singular
    /// value decomposition factors but also apply to jacobians as well as
    /// radiance values.
    //-----------------------------------------------------------------------
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

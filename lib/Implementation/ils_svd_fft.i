%include "fp_common.i"
%{
#include "ils_svd_fft.h"
#include "sub_state_vector_array.h"
%}

%base_import(ils)

%import "dispersion.i"
%import "state_vector.i"
%import "array_ad.i"

%fp_shared_ptr(FullPhysics::IlsSvdFft);

namespace FullPhysics {

class IlsSvdFft : public Ils {
public:
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
    virtual boost::shared_ptr<Ils> clone() const;
};

}

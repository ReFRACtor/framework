%include "fp_common.i"
%{
#include "ils_fast_apply.h"
#include "sub_state_vector_array.h"
%}

%base_import(ils)

%import "sample_grid.i"
%import "state_vector.i"
%import "array_ad.i"

%fp_shared_ptr(FullPhysics::IlsFastApply);

namespace FullPhysics {

class IlsFastApply : public Ils {
public:
    IlsFastApply(const blitz::Array<double, 2>& Left_matrix_truncated,
                 const blitz::Array<double, 2>& Right_matrix_fourier_transforms_real,
                 const blitz::Array<double, 2>& Right_matrix_fourier_transforms_imag,
                 const blitz::Array<int, 1>& Center_freq_indices,
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
    virtual SpectralDomain pixel_grid() const;
    virtual DoubleWithUnit high_res_extension() const;
    virtual void high_res_extension(const DoubleWithUnit& extension);
    virtual boost::shared_ptr<Ils> clone() const;
};

}

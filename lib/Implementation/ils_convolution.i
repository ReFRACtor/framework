%include "common.i"

%{
#include "ils_convolution.h"
#include "sub_state_vector_array.h"
%}

%base_import(ils_imp_base)
%import "ils_function.i"

%fp_shared_ptr(FullPhysics::IlsConvolution);

namespace FullPhysics {
class IlsConvolution : public IlsImpBase {
public:
    IlsConvolution(const boost::shared_ptr<SampleGrid>& Disp,
                   const boost::shared_ptr<IlsFunction>& Ils_func,
                   const DoubleWithUnit& Ils_half_width = DoubleWithUnit(20, units::inv_cm));
    IlsConvolution(const boost::shared_ptr<SampleGrid>& Disp,
                   const boost::shared_ptr<IlsFunction>& Ils_func,
                   double Ils_half_width);
    virtual blitz::Array<double, 1> apply_ils
    (const blitz::Array<double, 1>& High_resolution_wave_number,
     const blitz::Array<double, 1>& High_resolution_radiance,
     const std::vector<int>& Pixel_list) const;
    virtual ArrayAd<double, 1> apply_ils
    (const blitz::Array<double, 1>& High_resolution_wave_number,
     const ArrayAd<double, 1>& High_resolution_radiance,
     const std::vector<int>& Pixel_list) const;
    virtual boost::shared_ptr<Ils> clone() const;
    %python_attribute(ils_function, boost::shared_ptr<IlsFunction>);
};
}

%include "common.i"

%{
#include "mapping_scale.h"
%}

%base_import(mapping_imp_base)
%import "array_ad.i"
%fp_shared_ptr(FullPhysics::MappingScale);


namespace FullPhysics {
class MappingScale : public MappingImpBase  {
public:
    const ArrayAd<double, 1> fm_view(ArrayAd<double, 1> const& updated_coeff) const;
    const ArrayAd<double, 1> retrieval_init(ArrayAd<double, 1> const& initial_coeff) const;
    std::string name() const;
};
}
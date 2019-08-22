%include "fp_common.i"

%{
#include "mapping_offset.h"
%}

%base_import(mapping_imp_base)
%import "array_ad.i"
%import "pressure.i"
%fp_shared_ptr(FullPhysics::MappingOffset);


namespace FullPhysics {
class MappingOffset : public MappingImpBase {
public:
    const ArrayAd<double, 1> fm_view(ArrayAd<double, 1> const& updated_coeff,
        const boost::shared_ptr<Pressure>& updated_press) const;
    const ArrayAd<double, 1> retrieval_init(ArrayAd<double, 1> const& initial_coeff) const;
    std::string name() const;
};
}

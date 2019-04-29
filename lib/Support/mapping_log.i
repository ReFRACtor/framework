%include "common.i"

%{
#include "mapping_log.h"
%}

%base_import(mapping)
%import "array_ad.i"
%fp_shared_ptr(FullPhysics::MappingLog);


namespace FullPhysics {
class MappingLog : public Mapping {
public:
    const ArrayAd<double, 1> fm_view(ArrayAd<double, 1> const& coeff) const;
    const ArrayAd<double, 1> retrieval_view(ArrayAd<double, 1> const& coeff) const;
    std::string name() const;
};
}
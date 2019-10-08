%include "fp_common.i"

%{
#include "mapping_scale.h"
%}

%base_import(mapping)
%import "array_ad.i"
%import "pressure.i"
%fp_shared_ptr(FullPhysics::MappingScale);


namespace FullPhysics {
class MappingScale : public Mapping {
public:
    MappingScale(double Scale, blitz::Array<double, 1> Scalee);
    virtual ~MappingScale() {};
    
    virtual const ArrayAd<double, 1> fm_view(ArrayAd<double, 1> const& updated_coeff) const;
    virtual const ArrayAd<double, 1> retrieval_init(ArrayAd<double, 1> const& initial_coeff) const;
    virtual std::string name() const;
    virtual boost::shared_ptr<Mapping> clone() const;
};
}

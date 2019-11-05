%include "fp_common.i"

%{
#include "mapping_linear.h"
%}

%base_import(mapping)
%import "array_ad.i"
%import "pressure.i"
%fp_shared_ptr(FullPhysics::MappingLinear);


namespace FullPhysics {

%feature("notabstract") MappingLinear;

class MappingLinear : public Mapping  {
public:
    MappingLinear();
    virtual ~MappingLinear() {};
    virtual const ArrayAd<double, 1> fm_view(ArrayAd<double, 1> const& updated_coeff) const;
    virtual const ArrayAd<double, 1> retrieval_init(ArrayAd<double, 1> const& initial_coeff) const;
    virtual std::string name() const;
    virtual boost::shared_ptr<Mapping> clone() const;
};
}

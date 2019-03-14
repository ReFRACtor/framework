%include "common.i"

%{
#include "mapping.h"
%}

%base_import(generic_object)
%import "array_ad.i"
%fp_shared_ptr(FullPhysics::Mapping);


namespace FullPhysics {
class Mapping {
public:
    const ArrayAd<double, 1>& apply(ArrayAd<double, 1> const& coeff) const;
    AutoDerivative<double> invert(AutoDerivative<double> coeff_i) const;
};
}
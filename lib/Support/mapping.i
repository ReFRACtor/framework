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
    const ArrayAd<double, 1> apply(ArrayAd<double, 1> const& coeff) const;
    const blitz::Array<double, 1> apply(blitz::Array<double, 1> const& coeff) const;
    AutoDerivative<double> apply_element(AutoDerivative<double> coeff_i) const;
    const ArrayAd<double, 1> invert(ArrayAd<double, 1> const& coeff) const;
    const blitz::Array<double, 1> invert(blitz::Array<double, 1> const& coeff) const;
    AutoDerivative<double> invert_element(AutoDerivative<double> coeff_i) const;
    std::string name() const;
};
}
%include "fp_common.i"

%{
#include "rayleigh.h"
%}

%base_import(generic_object)

%import "array_ad.i"

%fp_shared_ptr(FullPhysics::Rayleigh);

namespace FullPhysics {
class Rayleigh: public GenericObject {
public:
    virtual ArrayAd<double, 1> optical_depth_each_layer(double wn, int spec_index) const = 0;

    virtual boost::shared_ptr<Rayleigh> clone() const;

    virtual void print(std::ostream& Os) const;
    std::string print_to_string() const;
    std::string print_parent() const;
    %pickle_serialization();
};
}

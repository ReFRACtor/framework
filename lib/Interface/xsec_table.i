%include "fp_common.i"

%{
#include "xsec_table.h"
%}

%base_import(generic_object)

%import "array_ad.i"
%import "double_with_unit.i"

%fp_shared_ptr(FullPhysics::XSecTable);

// Allow these classes to be derived from in Python.
%feature("director") TemperatureImpBase;

namespace FullPhysics {
class XSecTable: public GenericObject {
public:
    virtual ArrayAd<double, 2> optical_depth_each_layer(DoubleWithUnit grid_point, ArrayAd<double, 1> gas_density_levels, ArrayAd<double, 1> temperature_levels) const = 0;

    virtual boost::shared_ptr<XSecTable> clone() const = 0;

    virtual void print(std::ostream& Os) const;
    std::string print_to_string() const;
    %pickle_serialization();
};
}

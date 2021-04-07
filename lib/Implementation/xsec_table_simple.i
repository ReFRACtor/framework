%include "fp_common.i"

%{
#include "xsec_table_simple.h"
%}

%base_import(xsec_table_imp_base)

%fp_shared_ptr(FullPhysics::XSecTableSimple);

namespace FullPhysics {

class XSecTableSimple: public XSecTableImpBase {
public:
    XSecTableSimple(const ArrayWithUnit<double, 1>& Spectral_grid, const ArrayWithUnit<double, 2>& XSec_values, double Conversion_factor);

    virtual ArrayAdWithUnit<double, 1> optical_depth_each_layer_unweighted(DoubleWithUnit spectral_point, ArrayAdWithUnit<double, 1> gas_density_levels, ArrayAdWithUnit<double, 1> temperature_levels) const;

    virtual boost::shared_ptr<XSecTable> clone() const;

    virtual void print(std::ostream& Os) const;
    %pickle_serialization();
};

}

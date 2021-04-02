%include "fp_common.i"

%{
#include "xsec_table_temp_dep.h"
%}

%base_import(xsec_table_imp_base)

%fp_shared_ptr(FullPhysics::XSecTableTempDep);

namespace FullPhysics {

class XSecTableTempDep: public XSecTableImpBase {
public:
    XSecTableTempDep(const ArrayWithUnit<double, 1>& Spectral_grid, const blitz::Array<double, 2>& XSec_values, double Conversion_factor);

    virtual ArrayAd<double, 1> optical_depth_each_layer_unweighted(DoubleWithUnit spectral_point, ArrayAd<double, 1> gas_density_levels, ArrayAd<double, 1> temperature_levels) const;

    virtual boost::shared_ptr<XSecTable> clone() const;

    virtual void print(std::ostream& Os) const;
    %pickle_serialization();
};

}

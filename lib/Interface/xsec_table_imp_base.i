%include "fp_common.i"

%{
#include "xsec_table_imp_base.h"
%}

%base_import(xsec_table)

%import "array_with_unit.i"

%fp_shared_ptr(FullPhysics::XSecTableImpBase);

namespace FullPhysics {

// Allow these classes to be derived from in Python.
%feature("director") XSecTableImpBase;

class XSecTableImpBase: public XSecTable {
public:
    XSecTableImpBase(const ArrayWithUnit<double, 1>& Spectral_grid, const ArrayWithUnit<double, 2>& XSec_values, double Conversion_factor);
    virtual const ArrayWithUnit<double, 1> spectral_grid() const;
    virtual const DoubleWithUnit cross_section_value(DoubleWithUnit& spectral_point) const;
    virtual const ArrayWithUnit<double, 1> cross_section_coefficients(DoubleWithUnit& spectral_point) const;

    // From XSecTable, included for director usage
    virtual ArrayAdWithUnit<double, 1> optical_depth_each_layer_unweighted(DoubleWithUnit spectral_point, ArrayAdWithUnit<double, 1> gas_density_levels, ArrayAdWithUnit<double, 1> temperature_levels) const = 0;
    virtual boost::shared_ptr<XSecTable> clone() const = 0;

    virtual void print(std::ostream& Os) const;
    %pickle_serialization();
};

}
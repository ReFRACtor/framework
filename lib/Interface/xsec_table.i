%include "fp_common.i"

%{
#include "xsec_table.h"
%}

%base_import(generic_object)

%import "array_ad_with_unit.i"
%import "double_with_unit.i"

%fp_shared_ptr(FullPhysics::XSecTable);

namespace FullPhysics {

// Allow these classes to be derived from in Python.
%feature("director") XSecTable;

class XSecTable: public GenericObject {
public:
  XSecTable();
  virtual ArrayAdWithUnit<double, 1> optical_depth_each_layer_unweighted(DoubleWithUnit spectral_point, ArrayAdWithUnit<double, 1> gas_density_levels, ArrayAdWithUnit<double, 1> temperature_levels) const = 0;

  virtual boost::shared_ptr<XSecTable> clone() const = 0;

  virtual std::string desc() const;
  std::string print_to_string() const;
  std::string print_parent() const;
  %pickle_serialization();
};
}

%template(vector_xsec_table) std::vector<boost::shared_ptr<FullPhysics::XSecTable> >;

// Extra code for handling boost serialization/python pickle of
// director classes
%fp_director_serialization(xsec_table, XSecTable)

// List of things "import *" will include
%python_export("XSecTable", "vector_xsec_table");

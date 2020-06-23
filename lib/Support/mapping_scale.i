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
  virtual boost::shared_ptr<Mapping> clone();
  %python_attribute(initial_scale_factor, double);
  %python_attribute(scalee, blitz::Array<double, 1>);
  %pickle_serialization();
};
}

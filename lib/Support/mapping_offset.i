%include "fp_common.i"

%{
#include "mapping_offset.h"
%}

%base_import(mapping)
%import "array_ad.i"
%import "pressure.i"
%fp_shared_ptr(FullPhysics::MappingOffset);


namespace FullPhysics {

class MappingOffset : public Mapping {
public:
  MappingOffset(double Offset, blitz::Array<double, 1> Offsetee);
  virtual boost::shared_ptr<Mapping> clone();
  %python_attribute(initial_offset, double);
  %python_attribute(offsetee, blitz::Array<double, 1>);
  %pickle_serialization();
};
}

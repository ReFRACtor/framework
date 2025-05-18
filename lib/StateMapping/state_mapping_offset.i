%include "fp_common.i"

%{
#include "state_mapping_offset.h"
%}

%base_import(state_mapping)
%import "array_ad.i"
%import "pressure.i"
%fp_shared_ptr(FullPhysics::StateMappingOffset);


namespace FullPhysics {

class StateMappingOffset : public StateMapping {
public:
  StateMappingOffset(double Offset, blitz::Array<double, 1> Offsetee);
  virtual boost::shared_ptr<StateMapping> clone() const;
  %python_attribute(initial_offset, double);
  %python_attribute(offsetee, blitz::Array<double, 1>);
  %pickle_serialization();
};
}

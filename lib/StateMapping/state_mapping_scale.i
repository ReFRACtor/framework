%include "fp_common.i"

%{
#include "state_mapping_scale.h"
%}

%base_import(state_mapping)
%import "array_ad.i"
%import "pressure.i"
%fp_shared_ptr(FullPhysics::StateMappingScale);

namespace FullPhysics {
class StateMappingScale : public StateMapping {
public:
  StateMappingScale(double Scale, blitz::Array<double, 1> Scalee);
  virtual boost::shared_ptr<StateMapping> clone() const;
  %python_attribute(initial_scale_factor, double);
  %python_attribute(scalee, blitz::Array<double, 1>);
  %pickle_serialization();
};
}

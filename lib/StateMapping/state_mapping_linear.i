%include "fp_common.i"

%{
#include "state_mapping_linear.h"
%}

%base_import(state_mapping)
%import "array_ad.i"
%import "pressure.i"
%fp_shared_ptr(FullPhysics::StateMappingLinear);

namespace FullPhysics {

class StateMappingLinear : public StateMapping  {
public:
  StateMappingLinear();
  virtual boost::shared_ptr<StateMapping> clone() const;
  %pickle_serialization();
};
}

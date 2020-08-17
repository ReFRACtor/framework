%include "fp_common.i"

%{
#include "state_mapping_at_indexes.h"
%}

%base_import(state_mapping)
%import "array_ad.i"
%import "pressure.i"
%fp_shared_ptr(FullPhysics::StateMappingAtIndexes);

namespace FullPhysics {

class StateMappingAtIndexes : public StateMapping  {
public:
  StateMappingAtIndexes(const blitz::Array<int, 1>& indexes);
  StateMappingAtIndexes(const blitz::Array<bool, 1>& flags);
  virtual boost::shared_ptr<StateMapping> clone();
  %pickle_serialization();
};
}

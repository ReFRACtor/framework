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
  StateMappingAtIndexes();
  virtual boost::shared_ptr<StateMapping> clone();
  %pickle_serialization();
};
}

%include "fp_common.i"

%{
#include "state_mapping_log.h"
%}

%base_import(state_mapping)
%import "array_ad.i"
%import "pressure.i"
%fp_shared_ptr(FullPhysics::StateMappingLog);

namespace FullPhysics {
class StateMappingLog : public StateMapping {
public:
  StateMappingLog();
  virtual boost::shared_ptr<StateMapping> clone();
  %pickle_serialization();
};
}

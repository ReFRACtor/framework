%include "fp_common.i"

%{
#include "mapping_log.h"
%}

%base_import(mapping)
%import "array_ad.i"
%import "pressure.i"
%fp_shared_ptr(FullPhysics::MappingLog);

namespace FullPhysics {
class MappingLog : public Mapping {
public:
  MappingLog();
  virtual boost::shared_ptr<Mapping> clone();
  %pickle_serialization();
};
}

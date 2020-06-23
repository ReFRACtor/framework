// -*- mode: c++; -*-
// (Not really c++, but closest emacs mode)

%include "fp_common.i"

%{
#include "mapping_linear.h"
%}

%base_import(mapping)
%import "array_ad.i"
%import "pressure.i"
%fp_shared_ptr(FullPhysics::MappingLinear);

namespace FullPhysics {

class MappingLinear : public Mapping  {
public:
  MappingLinear();
  virtual boost::shared_ptr<Mapping> clone();
  %pickle_serialization();
};
}

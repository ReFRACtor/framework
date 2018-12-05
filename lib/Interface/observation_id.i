// -*- mode: c++; -*-
// (Not really c++, but closest emacs mode)

%include "common.i"
%{
#include "observation_id.h"
%}

%shared_ptr(FullPhysics::ObservationId)
namespace FullPhysics {
class ObservationId {
public:
  virtual ~ObservationId();
};
}




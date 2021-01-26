%include "fp_common.i"

%{
#include "state_vector_observer.h"
%}

%base_import(state_vector)
%base_import(observer)
    
%fp_shared_ptr(FullPhysics::StateVectorObserver);

%feature("director") FullPhysics::StateVectorObserver;

namespace FullPhysics {

class StateVectorObserver : public Observer<StateVector> {
public:
  virtual void state_vector_name(const StateVector& Sv, 
                                 blitz::Array<std::string, 1>& Sv_name) const;
  %python_attribute_with_set(sv_name, std::vector<std::string>);
  %pickle_serialization();
};

}

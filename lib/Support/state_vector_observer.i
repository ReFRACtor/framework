%include "fp_common.i"

%{
#include "state_vector_observer.h"
%}

%base_import(state_vector)
%base_import(observer)
    
%fp_shared_ptr(FullPhysics::StateVectorObserver);

namespace FullPhysics {

%nodefaultctor StateVectorObserver;

class StateVectorObserver : public Observer<StateVector> {
public:
  virtual ~StateVectorObserver();
  std::string print_to_string() const;
  virtual void mark_used(const StateVector& Sv, 
			 blitz::Array<bool, 1>& Used) const;
  virtual void state_vector_name(const StateVector& Sv, 
			  blitz::Array<std::string, 1>& Sv_name) const;
  %pickle_serialization();
};

}

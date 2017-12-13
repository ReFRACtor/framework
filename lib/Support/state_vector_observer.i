%include "common.i"

%{
#include "sub_state_vector_observer.h"
%}

%base_import(state_vector)
    
%fp_shared_ptr(FullPhysics::StateVectorObserver);
%fp_shared_ptr(FullPhysics::Observer<FullPhysics::StateVector>);

// Do this so we can derive from this and have it able to be used by the C++ code
// Defined here since rename does not like being inside of a namespace
%feature("director") FullPhysics::Observer<FullPhysics::StateVector>;
%rename(ObserverStateVector) FullPhysics::Observer<FullPhysics::StateVector>;

namespace FullPhysics {

class FullPhysics::Observer<FullPhysics::StateVector> {
public:
  virtual ~Observer<FullPhysics::StateVector>();
  virtual void notify_add(StateVector& Obs);
  virtual void notify_remove(StateVector& Obs);
  virtual void notify_update(const StateVector& Obs);
};

%nodefaultctor StateVectorObserver;

class StateVectorObserver : public Observer<StateVector> {
public:
  virtual ~StateVectorObserver();
  std::string print_to_string() const;
  virtual void mark_used(const StateVector& Sv, 
			 blitz::Array<bool, 1>& Used) const;
  virtual void state_vector_name(const StateVector& Sv, 
			  blitz::Array<std::string, 1>& Sv_name) const;
};

}

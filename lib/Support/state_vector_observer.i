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
  virtual void notify_update(const StateVector& Observed_object);
  virtual void notify_add(StateVector& Observed_object);
  virtual void notify_remove(StateVector& Observed_object);
  %python_attribute_with_set(sv_name, std::vector<std::string>);
  %pickle_serialization();
};
}

// Extra code for handling boost serialization/python pickle of
// director classes
%fp_director_serialization(state_vector_observer, StateVectorObserver);

// List of things "import *" will include
%python_export("StateVectorObserver");


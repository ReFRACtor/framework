// -*- mode: c++; -*-
// (Not really c++, but closest emacs mode)

%include "fp_common.i"

%{
#include "state_vector.h"
%}

%base_import(generic_object)
%import "observer.i"
%import "array_ad.i"

%fp_shared_ptr(FullPhysics::StateVector);
%fp_shared_ptr(FullPhysics::Observable<FullPhysics::StateVector>);
%fp_shared_ptr(FullPhysics::Observer<FullPhysics::StateVector>);

%feature("director") FullPhysics::Observer<FullPhysics::StateVector>;

namespace FullPhysics {
%template(ObserverStateVector) FullPhysics::Observer<FullPhysics::StateVector>;
%template(ObservableStateVector) FullPhysics::Observable<FullPhysics::StateVector>;

class StateVector: public Observable<StateVector> {
public:
  StateVector();
    
  virtual void add_observer(Observer<StateVector>& Obs);
  virtual void remove_observer(Observer<StateVector>& Obs);
  virtual void remove_observer(Observer<StateVector>& Obs, bool Keep_state_when_removed);
  void clear_observers(bool Keep_state_when_removed=false);
  %python_attribute(keep_state_when_removed, bool);

  %python_attribute(state, blitz::Array<double, 1>);
  %python_attribute(state_with_derivative, ArrayAd<double, 1>);
  %python_attribute(state_covariance, blitz::Array<double, 2>);
  

  void update_state(const blitz::Array<double, 1>& X);
  void update_state(const blitz::Array<double, 1>& X,
		    const blitz::Array<double, 2>& Cov);
  void update_state(const ArrayAd<double, 1>& X,
		    const blitz::Array<double, 2>& Cov);
    
  %python_attribute_with_set(observer_claimed_size, int);

  %python_attribute_with_set(default_state_vector_name, std::vector<std::string>);

  std::string print_to_string() const;

  %extend {
    std::vector<std::string> _state_vector_name() const 
    {
      blitz::Array<std::string, 1> sn($self->state_vector_name());
      std::vector<std::string> res;
      for(int i = 0; i < sn.extent(blitz::firstDim); ++i) {
	res.push_back(sn(i));
      }
      return res;
    }
  }

    %pythoncode {
        @property
        def state_vector_name(self):
            return self._state_vector_name()
    }
  %pickle_serialization();
};
}

// This silences  warning messages. Not sure if the message are important
// or not, but good to have them be quiet
%rename(ObserverStateVector2) FullPhysics::Observer<FullPhysics::StateVector>;

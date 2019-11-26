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
    virtual ~StateVector();
    StateVector();
    
    virtual void add_observer(Observer<StateVector>& Obs);
    virtual void remove_observer(Observer<StateVector>& Obs);

    %python_attribute(state, blitz::Array<double, 1>);
    %python_attribute(state_with_derivative, ArrayAd<double, 1>);
    %python_attribute(state_covariance, blitz::Array<double, 2>);

    void update_state(const blitz::Array<double, 1>& X);
    void update_state(const blitz::Array<double, 1>& X, const blitz::Array<double, 2>& Cov);
    
    %python_attribute(used_flag, blitz::Array<bool, 1>);
    %python_attribute_with_set(observer_claimed_size, int);

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

        // Support shared pointers from other languages for add/remove observer
        void add_observer(const boost::shared_ptr<Observer<StateVector> >& Obs) {
            $self->add_observer(*Obs);
        }

        void remove_observer(const boost::shared_ptr<Observer<StateVector> >& Obs) {
            $self->remove_observer(*Obs);
        }

    }

    %pythoncode {
        @property
        def state_vector_name(self):
            return self._state_vector_name()
    }

};
}

// This silences  warning messages. Not sure if the message are important
// or not, but good to have them be quiet
%rename(ObserverStateVector2) FullPhysics::Observer<FullPhysics::StateVector>;

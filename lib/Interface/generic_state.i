%include "fp_common.i"

%{
#include "generic_state.h"
%}

%base_import(observer)
%base_import(state_vector_observer)
%base_import(sub_state_vector_array)
%base_import(state_vector)
%import "auto_derivative.i"
%import "sub_state_vector_array.i"

%fp_shared_ptr(FullPhysics::GenericState)

%fp_shared_ptr(FullPhysics::Observable<FullPhysics::GenericState>)
%fp_shared_ptr(FullPhysics::Observer<FullPhysics::GenericState>)

%template(ObservableGenericState) FullPhysics::Observable<FullPhysics::GenericState>;
%template(ObserverGenericState) FullPhysics::Observer<FullPhysics::GenericState>;

namespace FullPhysics {
class GenericState : virtual public StateVectorObserver, public Observable<GenericState> {
public:
  virtual void add_observer(Observer<GenericState>& Obs);
  virtual void remove_observer(Observer<GenericState>& Obs);
  virtual boost::shared_ptr<GenericState> clone() const = 0;
  virtual std::string desc() const;
  std::string print_to_string() const;
  std::string print_parent() const;
  %pickle_serialization();
};
}

%fp_shared_ptr(FullPhysics::GenericStateImpBase);
%fp_shared_ptr(FullPhysics::SubStateVectorArray<FullPhysics::GenericState>)

namespace FullPhysics {
%template(SubStateVectorArrayGenericState) FullPhysics::SubStateVectorArray<FullPhysics::GenericState>;
  
%feature("director") GenericStateImpBase;
  
class GenericStateImpBase : public SubStateVectorArray<GenericState> {
public:
  GenericStateImpBase();
  virtual boost::shared_ptr<GenericState> clone() const;
  %sub_state_virtual_func(GenericState);
  virtual std::string desc() const;
  std::string print_to_string() const;
  std::string print_parent() const;
  %pickle_serialization();
};
 
}

// Extra code for handling boost serialization/python pickle of
// director classes
%fp_director_serialization(generic_state, GenericStateImpBase)

// List of things "import *" will include
%python_export("GenericState", "GenericStateImpBase");


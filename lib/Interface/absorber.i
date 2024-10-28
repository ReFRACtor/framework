%include "fp_common.i"
%{
#include "sub_state_vector_array.h"
#include "absorber.h"
#include "temperature.h"
#include "altitude.h"
%}

%base_import(state_vector_observer)
%base_import(sub_state_vector_observer)
%base_import(observer)

%import "sub_state_vector_array.i"
%import "sub_state_vector_observer.i"
%import "pressure.i"
%import "temperature.i"
%import "altitude.i"
%import "absorber_vmr.i"

namespace FullPhysics {
  class Absorber;
}

%nodefaultctor FullPhysics::SubStateVectorArray<FullPhysics::Absorber>;
%fp_shared_ptr(FullPhysics::Absorber)
%fp_shared_ptr(FullPhysics::Observable<FullPhysics::Absorber>)
%fp_shared_ptr(FullPhysics::Observer<FullPhysics::Absorber>)
%fp_shared_ptr(FullPhysics::SubStateVectorArray<FullPhysics::Absorber>);

namespace FullPhysics {
%template(ObservableAbsorber) FullPhysics::Observable<Absorber>;
%template(ObserverAbsorber) FullPhysics::Observer<Absorber>;

// Allow these classes to be derived from in Python.
%feature("director") Absorber;

// Note, a class that is derived from in python needs to declare every virtual function that
// can be called on it, even if all that happens is the base class
// to a director is called. This is because this class is used to
// create the SwigDirector class, and this class needs each of the member functions to
// direct things properly. It is *not* necessary to add these function to the underlying
// C++, only that you declare them here.

class Absorber : virtual public StateVectorObserver, 
		 public Observable<Absorber> {
public:
  virtual ~Absorber();
  virtual std::string desc() const;
  virtual void add_observer(Observer<Absorber>& Obs);
  virtual void remove_observer(Observer<Absorber>& Obs);
  std::string print_to_string() const;
  std::string print_parent() const;
  %python_attribute(number_species, virtual int);
  virtual std::string gas_name(int Species_index) const = 0;
  virtual int gas_index(const std::string& Name) const;

  virtual ArrayAdWithUnit<double, 1> total_air_number_density_layer(int spec_index) const = 0;
  virtual ArrayAdWithUnit<double, 2> gas_number_density_layer(int spec_index) const = 0;
 
  virtual ArrayAd<double, 2> optical_depth_each_layer(double wn, int spec_index) const = 0;
  virtual boost::shared_ptr<AbsorberVmr> absorber_vmr(const std::string& gas_name) const = 0;
  virtual boost::shared_ptr<Absorber> clone() const = 0;

  // Functions so StateVectorObserver base class, see note above
  virtual void notify_update(const StateVector& Sv);
  virtual void state_vector_name(const StateVector& Sv, 
			  blitz::Array<std::string, 1>& Sv_name) const;
  virtual void notify_add(StateVector& Observed_object);
  virtual void notify_remove(StateVector& Observed_object);
  %python_attribute(subobject_list, std::vector<boost::shared_ptr<GenericObject> >);
  %pickle_serialization();
};
}

// Extra code for handling boost serialization/python pickle of
// director classes
%fp_director_serialization(absorber, Absorber)

// List of things "import *" will include
%python_export("Absorber");

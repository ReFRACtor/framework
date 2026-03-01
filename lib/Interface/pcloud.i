%include "fp_common.i"

%{
#include "pcloud.h"
%}

%base_import(state_vector_observer)
%import "array_with_unit.i"
%import "auto_derivative_with_unit.i"

%fp_shared_ptr(FullPhysics::Pcloud)

%fp_shared_ptr(FullPhysics::Observable<FullPhysics::Pcloud>)
%fp_shared_ptr(FullPhysics::Observer<FullPhysics::Pcloud>)

%template(ObservablePcloud) FullPhysics::Observable<FullPhysics::Pcloud>;
%template(ObserverPcloud) FullPhysics::Observer<FullPhysics::Pcloud>;

namespace FullPhysics {
class Pcloud : virtual public StateVectorObserver, public Observable<Pcloud> {
public:
  virtual void add_observer(Observer<Pcloud>& Obs);
  virtual void remove_observer(Observer<Pcloud>& Obs);
  virtual AutoDerivativeWithUnit<double> pressure_cloud(int sensor_index) const = 0;
  virtual boost::shared_ptr<Pcloud> clone() const = 0;
  virtual std::string desc() const;
  std::string print_to_string() const;
  std::string print_parent() const;
  %pickle_serialization();
};
}

// List of things "import *" will include
%python_export("Pcloud");


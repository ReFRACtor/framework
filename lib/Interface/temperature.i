// -*- mode: c++; -*-
// (Not really c++, but closest emacs mode)

%include "fp_common.i"

%{
#include "temperature.h"
%}

%base_import(state_vector_observer)
%import "array_with_unit.i"
%import "auto_derivative_with_unit.i"
%import "pressure.i"

%fp_shared_ptr(FullPhysics::Temperature)
namespace FullPhysics {
  class Temperature;
}

%fp_shared_ptr(FullPhysics::Observable<FullPhysics::Temperature>)
%fp_shared_ptr(FullPhysics::Observer<FullPhysics::Temperature>)

namespace FullPhysics {
%template(ObservableTemperature) FullPhysics::Observable<Temperature>;
%template(ObserverTemperature) FullPhysics::Observer<Temperature>;

class Temperature : virtual public StateVectorObserver, 
		    public Observable<Temperature> {
public:
  virtual ~Temperature();
  virtual void add_observer(Observer<Temperature>& Obs); 
  virtual void remove_observer(Observer<Temperature>& Obs);
  %python_attribute(important_pressure_level, virtual ArrayWithUnit<double, 1>);
  virtual AutoDerivativeWithUnit<double> 
  temperature(const AutoDerivativeWithUnit<double>& Press) const = 0;
  virtual ArrayAdWithUnit<double, 1> temperature_grid(const Pressure& P) const;
  virtual boost::shared_ptr<Temperature> clone() const = 0;
  std::string print_to_string() const;
  %pickle_serialization();
};
} 


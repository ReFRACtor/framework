// -*- mode: c++; -*-
// (Not really c++, but closest emacs mode)

%include "fp_common.i"

%{
#include "pressure.h"
%}

%base_import(state_vector_observer)
%import "array_ad_with_unit.i"
%import "auto_derivative_with_unit.i"

%fp_shared_ptr(FullPhysics::Pressure)
namespace FullPhysics {
  class Pressure;
}

%fp_shared_ptr(FullPhysics::Observable<FullPhysics::Pressure>)
%fp_shared_ptr(FullPhysics::Observer<FullPhysics::Pressure>)

%template(ObservablePressure) FullPhysics::Observable<FullPhysics::Pressure>;
%template(ObserverPressure) FullPhysics::Observer<FullPhysics::Pressure>;

namespace FullPhysics {
// Allow this class to be derived from in Python.
%feature("director") Pressure;


class Pressure : virtual public StateVectorObserver, 
		 public Observable<Pressure> {
public:
  enum TypePreference {PREFER_INCREASING_PRESSURE=0,
    PREFER_DECREASING_PRESSURE=1};
  enum PressureGridType { INCREASING_PRESSURE=0,
    DECREASING_PRESSURE=1, NATIVE_ORDER=2};
  virtual ~Pressure();
  virtual void add_observer(Observer<Pressure>& Obs);
  virtual void remove_observer(Observer<Pressure>& Obs);
  virtual ArrayAdWithUnit<double, 1>
  pressure_grid(PressureGridType Gtype = INCREASING_PRESSURE) const = 0;
  %python_attribute(surface_pressure, AutoDerivativeWithUnit<double>)
  %python_attribute(surface_pressure_value, double)
  %python_attribute_abstract(type_preference, TypePreference)
  %python_attribute(number_layer, int)
  %python_attribute(number_level, int)
  %python_attribute(max_number_level, virtual int)
  %pickle_serialization();
  virtual boost::shared_ptr<Pressure> clone() const = 0;
  std::string print_to_string() const;
  %pickle_serialization();
};
}


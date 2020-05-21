// -*- mode: c++; -*-
// (Not really c++, but closest emacs mode)

%include "fp_common.i"

%{
#include "scalar_with_unit.h"
%}

%base_import(generic_object)
%import "unit.i"

namespace FullPhysics {
template<class T> class ScalarWithUnit : public GenericObject {
public:
  std::string print_to_string() const;
  ScalarWithUnit();
  ScalarWithUnit(const ScalarWithUnit<T>& V);
  ScalarWithUnit(T V, const Unit& U);
  ScalarWithUnit(T V, const std::string& U);
  ScalarWithUnit(T V);
  ScalarWithUnit<T>& operator*=(const ScalarWithUnit<T>& D);
  ScalarWithUnit<T>& operator/=(const ScalarWithUnit<T>& D);
  ScalarWithUnit<T>& operator+=(const ScalarWithUnit<T>& D);
  ScalarWithUnit<T>& operator-=(const ScalarWithUnit<T>& D);
  ScalarWithUnit<T> convert(const Unit& R) const;
  ScalarWithUnit<T> convert(const std::string& R) const;
  ScalarWithUnit<T> convert_wave(const Unit& R) const;
  ScalarWithUnit<T> convert_wave(const std::string& R) const;
  // There is a cyclic dependency here. Go ahead and just leave this
  // out for now, we can come up with a solution if we need to.
  // ScalarWithUnit<T> convert_wave(const Unit& R, const SpectralDomain& Sd) const;
  %extend {
    T _value() const {return $self->value; }
    void _value_set(T V) { $self->value = V;}
    Unit _units() const {return $self->units; }
    void _units_set(const FullPhysics::Unit& U) {$self->units = U;}
    ScalarWithUnit<T> __mul__(const ScalarWithUnit<T>& Y) 
    { return *$self * Y; }
    // Python 2 division operator name
    ScalarWithUnit<T> __div__(const ScalarWithUnit<T>& Y) 
    { return *$self / Y; }
    // Python 3 division operator name
    ScalarWithUnit<T> __truediv__(const ScalarWithUnit<T>& Y) 
    { return *$self / Y; }
    ScalarWithUnit<T> __add__(const ScalarWithUnit<T>& Y) 
    { return *$self + Y; }
    ScalarWithUnit<T> __sub__(const ScalarWithUnit<T>& Y) 
    { return *$self - Y; }
  }
  %pythoncode {
@property
def value(self):
  return self._value()

@value.setter
def value(self, val):
  self._value_set(val)

@property
def units(self):
  return self._units()

@units.setter
def units(self,val):
    self._units_set(val)
  }
};
}
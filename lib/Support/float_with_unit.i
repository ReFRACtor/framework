// -*- mode: c++; -*-
// (Not really c++, but closest emacs mode)

%include "fp_common.i"

%{
#include "float_with_unit.h"
%}

%base_import(generic_object)
%import "unit.i"
%fp_shared_ptr(FullPhysics::FloatWithUnit)
namespace FullPhysics {
class FloatWithUnit : public GenericObject {
public:
  std::string print_to_string() const;
  FloatWithUnit();
  FloatWithUnit(float V, const Unit& U);
  FloatWithUnit(float V, const std::string& U);
  FloatWithUnit(float V);
  FloatWithUnit& operator*=(const FloatWithUnit& D);
  FloatWithUnit& operator/=(const FloatWithUnit& D);
  FloatWithUnit& operator+=(const FloatWithUnit& D);
  FloatWithUnit& operator-=(const FloatWithUnit& D);
  FloatWithUnit convert(const Unit& R) const;
  FloatWithUnit convert(const std::string& R) const;
  FloatWithUnit convert_wave(const Unit& R) const;
  FloatWithUnit convert_wave(const std::string& R) const;
  // There is a cyclic dependency here. Go ahead and just leave this
  // out for now, we can come up with a solution if we need to.
  // FloatWithUnit convert_wave(const Unit& R, const SpectralDomain& Sd) const;
  %extend {
    float _value() const {return $self->value; }
    void _value_set(float V) { $self->value = V;}
    Unit _units() const {return $self->units; }
    void _units_set(const FullPhysics::Unit& U) {$self->units = U;}
    FloatWithUnit __mul__(const FloatWithUnit& Y) 
    { return *$self * Y; }
    // Python 2 division operator name
    FloatWithUnit __div__(const FloatWithUnit& Y) 
    { return *$self / Y; }
    // Python 3 division operator name
    FloatWithUnit __truediv__(const FloatWithUnit& Y) 
    { return *$self / Y; }
    FloatWithUnit __add__(const FloatWithUnit& Y) 
    { return *$self + Y; }
    FloatWithUnit __sub__(const FloatWithUnit& Y) 
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


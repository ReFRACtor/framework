// -*- mode: c++; -*-
// (Not really c++, but closest emacs mode)

%include "fp_common.i"

%{
#include "array_ad_with_unit.h"
%}

%import "array_ad.i"
%import "unit.i"

%fp_shared_ptr(FullPhysics::ArrayAdWithUnit<double, 1>);
%fp_shared_ptr(FullPhysics::ArrayAdWithUnit<double, 2>);
%fp_shared_ptr(FullPhysics::ArrayAdWithUnit<double, 3>);
%fp_shared_ptr(FullPhysics::ArrayAdWithUnit<double, 4>);

%pythoncode %{
import numpy as np
from .auto_derivative_with_unit import AutoDerivativeWithUnitDouble
from .auto_derivative import AutoDerivativeDouble
from .array_ad import (ArrayAd_double_1, ArrayAd_double_2, ArrayAd_double_3)
from .unit import Unit
%}

namespace FullPhysics {
template<class T, int D> class ArrayAdWithUnit: public GenericObject
{
public:
  ArrayAdWithUnit();
  ArrayAdWithUnit(const ArrayAd<T, D>& V, const Unit& U);
  ArrayAdWithUnit(const ArrayAd<T, D>& V, const std::string& U);
  ArrayAdWithUnit(const ArrayAd<T, D>& V);
  ArrayAdWithUnit<T, D> convert(const Unit& R) const;
  ArrayAdWithUnit<T, D> convert(const std::string& R) const;
  %python_attribute(rows, int)
  %python_attribute(cols, int)
  %python_attribute(depth, int)
  %python_attribute(is_constant, bool)
  %python_attribute(number_variable, int)
  void reference(const ArrayAdWithUnit<T, D>& V);
  %extend {
    FullPhysics::ArrayAd<T, D> _value() const { return $self->value;}
    void _value_set(const FullPhysics::ArrayAd<T, D>& V) { $self->value = V;}
    FullPhysics::Unit _units() const {return $self->units;}
    void _units_set(const FullPhysics::Unit& U) {$self->units = U;}
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

def __getitem__(self, index):
    sel_vals = self._value()[index]

    if isinstance(sel_vals, AutoDerivativeDouble):
        return AutoDerivativeWithUnitDouble(sel_vals, self.units)
    elif(isinstance(sel_vals, ArrayAd_double_1)):
        return ArrayAdWithUnit_double_1(sel_vals, self.units)
    elif(isinstance(sel_vals, ArrayAd_double_2)):
        return ArrayAdWithUnit_double_2(sel_vals, self.units)
    elif(isinstance(sel_vals, ArrayAd_double_3)):
        return ArrayAdWithUnit_double_3(sel_vals, self.units)
    else:
        raise NonImplementedError("__getitem__ limited to extracting slices of up to 3 dimensions, for type AutoDerivativeDouble")
  
def __setitem__(self, index, val):
    raise NotImplementedError("Setting values not yet implemented as it would require multiple data copies")
  }
  std::string print_to_string() const;
  std::string print_parent() const;
  %pickle_serialization();
};
%template(ArrayAdWithUnit_double_1) FullPhysics::ArrayAdWithUnit<double, 1>;
%template(ArrayAdWithUnit_double_2) FullPhysics::ArrayAdWithUnit<double, 2>;
%template(ArrayAdWithUnit_double_3) FullPhysics::ArrayAdWithUnit<double, 3>;
%template(ArrayAdWithUnit_double_4) FullPhysics::ArrayAdWithUnit<double, 4>;
}





// -*- mode: c++; -*-
// (Not really c++, but closest emacs mode)

%include "fp_common.i"

%{
#include "pressure_imp_base.h"
%}

%base_import(observer)
%base_import(pressure)
%base_import(sub_state_vector_array)
%base_import(state_vector)

%fp_shared_ptr(FullPhysics::PressureImpBase);
%fp_shared_ptr(FullPhysics::PressureImpBaseCache);
%fp_shared_ptr(FullPhysics::SubStateVectorArray<FullPhysics::Pressure>)

namespace FullPhysics {

%template(SubStateVectorArrayPressure) FullPhysics::SubStateVectorArray<Pressure>;

class PressureImpBaseCache : public CacheInvalidatedObserver {
public:
  PressureImpBaseCache();
  ArrayAdWithUnit<double, 1> pgrid;
  %pickle_serialization();
};
// Allow these classes to be derived from in Python.
%feature("director") PressureImpBase;

// Note, at least for SWIG 2.0.4 a class that is derived from in python 
// needs to declare every virtual function that can be called on it, even 
// if all that happens is the base class to a director is called. This is 
// because this class is used to create the SwigDirector class, and this 
// class needs each of the member functions to direct things properly. It 
// is *not* necessary to add these function to the underlying
// C++, only that you declare them here.
//
// If you miss something, then you will get something like a recursion
// error in python when a virtual function is used that isn't explicitly
// listed here.
//
// This seems like a bug in 2.0.4, if SWIG needs all the member functions
// it should know to make them itself. So perhaps a future version of SWIG
// won't have this same constraint. But for now, this is required.

class PressureImpBase: public SubStateVectorArray<Pressure> {
public:
  // From PressureImpBase
  virtual ~PressureImpBase();
  virtual ArrayAdWithUnit<double, 1> pressure_grid(Pressure::PressureGridType Gtype = Pressure::INCREASING_PRESSURE) const;
  %python_attribute_derived(type_preference, Pressure::TypePreference);
  %sub_state_virtual_func(Pressure);
  %pickle_serialization();
protected:
  mutable PressureImpBaseCache cache;
  virtual void calc_pressure_grid() const = 0;
  PressureImpBase();
  PressureImpBase(const blitz::Array<double, 1>& Coeff,
		  Pressure::TypePreference Tpref =
		  Pressure::PREFER_INCREASING_PRESSURE);
};
}


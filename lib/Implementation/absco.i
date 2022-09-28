// -*- mode: c++; -*-
// (Not really c++, but closest emacs mode)
%include <std_vector.i>
%include "fp_common.i"

%{
#include "absco.h"
%}

%import "double_with_unit.i"
%import "auto_derivative_with_unit.i"

%base_import(gas_absorption)

%fp_shared_ptr(FullPhysics::Absco);

namespace FullPhysics {

// Allow these classes to be derived from in Python.
%feature("director") Absco;

class Absco : public GasAbsorption {
public:
  int number_broadener_vmr(int Broadener_index) const;
  %python_attribute(number_layer, int)
  %python_attribute(number_temperature, int)
  virtual blitz::Array<double, 1> broadener_vmr_grid(int Broadner_index) const = 0;
  %python_attribute_abstract(pressure_grid, blitz::Array<double, 1>)
  %python_attribute_abstract(temperature_grid, blitz::Array<double, 2>)
  virtual double table_scale(double wn) const = 0;

  virtual DoubleWithUnit absorption_cross_section(double Wn, 
     const DoubleWithUnit& Press, 
     const DoubleWithUnit& Temp,
     const ArrayWithUnit<double, 1>& Broadener_vmr) const;

  virtual AutoDerivativeWithUnit<double> absorption_cross_section(double wn, 
    const DoubleWithUnit& press, 
    const AutoDerivativeWithUnit<double>& temp,
    const ArrayAdWithUnit<double, 1>& Broadener_vmr) const;

  virtual bool is_float() const = 0;

  virtual blitz::Array<double, 3> read_double(double wn) const = 0;
  virtual blitz::Array<float, 3> read_float(double wn) const = 0;
  virtual blitz::Array<double, 4> read_double_2b(double wn) const;
  virtual blitz::Array<float, 4> read_float_2b(double wn) const;

  %pickle_serialization();
};
}

%template(vector_absco) std::vector<boost::shared_ptr<FullPhysics::Absco> >;

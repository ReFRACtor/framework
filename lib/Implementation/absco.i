// -*- mode: c++; -*-
// (Not really c++, but closest emacs mode)
%include <std_vector.i>
%include "fp_common.i"
%{
#include "absco.h"
%}
%base_import(gas_absorption)
%import "double_with_unit.i"
%import "auto_derivative_with_unit.i"
%fp_shared_ptr(FullPhysics::Absco);
namespace FullPhysics {
class Absco : public GasAbsorption {
public:
  int number_broadener_vmr(int Broadener_index) const;
  %python_attribute(number_layer, int)
  %python_attribute(number_temperature, int)
  blitz::Array<double, 1> broadener_vmr_grid(int Broadner_index) const;
  %python_attribute_abstract(pressure_grid,blitz::Array<double, 1>)
  %python_attribute_abstract(temperature_grid,blitz::Array<double, 2>)
  virtual double table_scale(double wn) const;
  virtual DoubleWithUnit absorption_cross_section(double Wn, 
     const DoubleWithUnit& Press, 
     const DoubleWithUnit& Temp,
     const ArrayWithUnit<double, 1>& Broadener_vmr) const;
  virtual AutoDerivativeWithUnit<double>
  absorption_cross_section(double wn, 
    const DoubleWithUnit& press, 
    const AutoDerivativeWithUnit<double>& temp,
    const ArrayAdWithUnit<double, 1>& Broadener_vmr) const;
  %extend {
    blitz::Array<double, 3> read_double(double wn) const 
    { return $self->read<double, 3>(wn); }
    blitz::Array<float, 3> read_float(double wn) const 
    { return $self->read<float, 3>(wn); }
  }
  %pickle_serialization();
};
}

%template(vector_absco) std::vector<boost::shared_ptr<FullPhysics::Absco> >;

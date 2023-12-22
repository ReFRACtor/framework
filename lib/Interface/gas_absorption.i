// -*- mode: c++; -*-
// (Not really c++, but closest emacs mode)

%include <std_vector.i>
%include "fp_common.i"

%{#include "gas_absorption.h"
#include "gas_absorption.h"
%}

%base_import(generic_object)

%import "double_with_unit.i"
%import "auto_derivative_with_unit.i"
%import "array_with_unit.i"
%import "array_ad_with_unit.i"

%fp_shared_ptr(FullPhysics::GasAbsorption);

namespace FullPhysics {

// Allow these classes to be derived from in Python.
%feature("director") GasAbsorption;

class GasAbsorption : public GenericObject {
public:
  virtual ~GasAbsorption();
  virtual std::string desc() const;
  std::string print_to_string() const;
  virtual bool have_data(double wn) const = 0;
  %python_attribute(number_broadener, virtual int);
  virtual std::string broadener_name(int Broadner_index) const = 0;
  virtual DoubleWithUnit absorption_cross_section(double Wn, 
     const DoubleWithUnit& Press, 
     const DoubleWithUnit& Temp,
     const ArrayWithUnit<double, 1>& Broadener_vmr) const = 0;
  virtual AutoDerivativeWithUnit<double>
  absorption_cross_section(double wn, 
    const DoubleWithUnit& press, 
    const AutoDerivativeWithUnit<double>& temp,
    const ArrayAdWithUnit<double, 1>& Broadener_vmr) const = 0;
  virtual void print(std::ostream& Os) const;
  %pickle_serialization();
};
}

%template(vector_gas_absorption) std::vector<boost::shared_ptr<FullPhysics::GasAbsorption> >;

// Extra code for handling boost serialization/python pickle of
// director classes
%{
// Needed by code below, can't easily figure these names out
// automatically so just include here
#include "gas_absorption_wrap.h"
%}
%fp_director_serialization(GasAbsorption)

// List of things "import *" will include
%python_export("GasAbsorption");

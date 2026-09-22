%include <std_vector.i>
%include "fp_common.i"

%{
#include "absco.h"
%}

%import "double_with_unit.i"
%import "auto_derivative_with_unit.i"
%import "array_ad.i"
%import "array_with_unit.i"
%import "array_ad_with_unit.i"

%base_import(gas_absorption)

// Note we need the template before its first use, so we have all
// the typemaps in place.
%template(vector_absco) std::vector<boost::shared_ptr<FullPhysics::Absco> >;

// We also need the shared ptr after the templates, so we override the
// output type maps for std::vector<boost::shared_ptr<T> >
%fp_shared_ptr(FullPhysics::Absco);
%fp_shared_ptr(FullPhysics::AbscoInterpolator);

namespace FullPhysics {

// Allow these classes to be derived from in Python.
%feature("director") Absco;

class Absco : public GasAbsorption {
public:
  int number_broadener_vmr(int Broadener_index) const;
  virtual std::string desc() const;
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

// Helper class that calculates the absorption cross section for a
// fixed set of Pressure, Temperature, and Broadener VMR. Needs to be
// declared after Absco since its constructor references it.
class AbscoInterpolator : public GenericObject {
public:
  AbscoInterpolator(const boost::shared_ptr<Absco>& A,
		    const ArrayWithUnit<double, 1>& Press,
		    const ArrayAdWithUnit<double, 1>& Temp,
		    const ArrayAdWithUnit<double, 2>& Broadener_vmr);
  std::string print_to_string() const;
  std::string print_parent() const;
  blitz::Array<double, 1> absorption_cross_section_noderiv(double wn) const;
  ArrayAd<double, 1> absorption_cross_section_deriv(double wn) const;
};
}

// Extra code for handling boost serialization/python pickle of
// director classes
%fp_director_serialization(absco, Absco)


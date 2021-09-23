// -*- mode: c++; -*-
// (Not really c++, but closest emacs mode)
%include "fp_common.i"
%{
#include "pressure_sigma.h"
%}

%base_import(pressure_imp_base)
%import "auto_derivative.i"
%fp_shared_ptr(FullPhysics::PressureSigma);

namespace FullPhysics {
class PressureSigma : public PressureImpBase {
public:
  PressureSigma(const blitz::Array<double, 1>& A,
                const blitz::Array<double, 1>& B,
                double Surface_pressure,
		Pressure::TypePreference Tpref = Pressure::PREFER_INCREASING_PRESSURE);

  PressureSigma(const blitz::Array<double, 1>& Pressure_grid,
                double Surface_pressure,
		Pressure::TypePreference Tpref = Pressure::PREFER_INCREASING_PRESSURE);
  %python_attribute(surface_pressure_uncertainty, double)
  %rename(_v_set_surface_pressure) set_surface_pressure;
  void set_surface_pressure(const AutoDerivative<double>& Surface_pressure);
  %pythoncode {
    @property
    def surface_pressure(self):
      return self._v_surface_pressure()
      
    @surface_pressure.setter
    def surface_pressure(self, value):
       self._v_set_surface_pressure(value)
  }
  void set_levels_from_grid(const blitz::Array<double, 1>& Pressure_grid);
  virtual boost::shared_ptr<Pressure> clone() const;
  %python_attribute(sub_state_identifier, std::string);
  virtual std::string state_vector_name_i(int i) const;
  %python_attribute(a, blitz::Array<double, 1>)
  %python_attribute(b, blitz::Array<double, 1>)
  %pickle_serialization();
protected:
  virtual void calc_pressure_grid() const;
};
}

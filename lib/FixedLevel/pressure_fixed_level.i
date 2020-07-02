// -*- mode: c++; -*-
// (Not really c++, but closest emacs mode)
%include "fp_common.i"
%{
#include "pressure_fixed_level.h"
%}

%base_import(pressure_imp_base)
%import "pressure_level_input.i"

%fp_shared_ptr(FullPhysics::PressureFixedLevel);

namespace FullPhysics {
  class PressureFixedLevel;

class PressureFixedLevel : public PressureImpBase {
public:
  PressureFixedLevel(bool Pressure_flag, 
		     const boost::shared_ptr<PressureLevelInput>& Press_level,
		     double Surface_pressure);
  %python_attribute(surface_pressure_uncertainty, double);
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
  %python_attribute(number_active_level, int)
  %python_attribute(number_active_layer, int)
  %python_attribute(max_number_level, int)
  %python_attribute(pressure_active_levels, blitz::Array<double, 1>)
  virtual boost::shared_ptr<Pressure> clone() const;
  %python_attribute(sub_state_identifier, std::string);
  virtual std::string state_vector_name_i(int i) const;
  %pickle_serialization();
protected:
  virtual void calc_pressure_grid() const;
};
}

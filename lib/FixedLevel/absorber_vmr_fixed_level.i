// -*- mode: c++; -*-
// (Not really c++, but closest emacs mode)
%include "fp_common.i"
%{
#include "absorber_vmr_fixed_level.h"
%}

%base_import(absorber_vmr_imp_base)
%import "pressure.i"
%import "pressure_level_input.i"
%fp_shared_ptr(FullPhysics::AbsorberVmrFixedLevel)
namespace FullPhysics {
class AbsorberVmrFixedLevel : public AbsorberVmrImpBase {
public:
  AbsorberVmrFixedLevel(const boost::shared_ptr<Pressure>& Press,
	const boost::shared_ptr<PressureLevelInput>& Press_level,	   
			const blitz::Array<bool, 1>& Flag, 
			const blitz::Array<double, 1>& Vmr,
			const std::string& Gas_name);
  virtual boost::shared_ptr<AbsorberVmr> clone() const;
  %python_attribute(sub_state_identifier, std::string);
  virtual std::string state_vector_name_i(int i) const;
  %python_attribute(volume_mixing_ratio_level, blitz::Array<double, 1>);
  %python_attribute(volume_mixing_ratio_active_level, blitz::Array<double, 1>);
protected:
  virtual void calc_vmr() const;
};
}


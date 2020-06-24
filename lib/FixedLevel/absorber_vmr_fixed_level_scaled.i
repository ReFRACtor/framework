// -*- mode: c++; -*-
// (Not really c++, but closest emacs mode)
%include "fp_common.i"
%{
#include "absorber_vmr_fixed_level_scaled.h"
%}
%base_import(absorber_vmr_imp_base)
%import "pressure.i"
%import "pressure_level_input.i"

%fp_shared_ptr(FullPhysics::AbsorberVmrFixedLevelScaled)
namespace FullPhysics {
class AbsorberVmrFixedLevelScaled : public AbsorberVmrImpBase {
public:
  AbsorberVmrFixedLevelScaled(const boost::shared_ptr<Pressure>& Press,
        const boost::shared_ptr<PressureLevelInput>& Press_level,          
        const blitz::Array<double, 1>& Vmr,
        bool Used_flag,
        double Scale,                         
        const std::string& Gas_name);
  %python_attribute(sub_state_identifier, std::string);
  virtual std::string state_vector_name_i(int i) const;
  virtual boost::shared_ptr<AbsorberVmr> clone() const;
  %python_attribute(scale_factor, double);
  %python_attribute(scale_uncertainty, double);
  %pickle_serialization()
protected:
  virtual void calc_vmr() const;
};
}


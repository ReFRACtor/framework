%include "fp_common.i"
%{
#include "absorber_vmr_level_scaled.h"
%}

%base_import(absorber_vmr_level)
%import "pressure.i"
%fp_shared_ptr(FullPhysics::AbsorberVmrLevelScaled)

namespace FullPhysics {
class AbsorberVmrLevelScaled : public AbsorberVmrLevel {
public:
  AbsorberVmrLevelScaled(const boost::shared_ptr<Pressure>& Press,
                         const blitz::Array<double, 1>& Vmr_profile,
                         double Scale,                         
                         bool Scale_flag,
                         const std::string& Gas_name);

  virtual boost::shared_ptr<AbsorberVmr> clone() const;

  %python_attribute(scale_factor, double)
  %python_attribute(vmr_profile, blitz::Array<double, 1>)
  %python_attribute(pressure_profile, blitz::Array<double, 1>)

};
}


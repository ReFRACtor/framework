%include "fp_common.i"

%{
#include "absorber_vmr_level.h"
%}

%base_import(absorber_vmr_imp_base)
%base_import(state_mapping)
%base_import(state_mapping_linear)
%import "pressure.i"

%fp_shared_ptr(FullPhysics::AbsorberVmrLevel)

namespace FullPhysics {
class AbsorberVmrLevel : public AbsorberVmrImpBase {
public:
  AbsorberVmrLevel(const boost::shared_ptr<Pressure>& Press,
                   const blitz::Array<double, 1>& Vmr,
                   const blitz::Array<bool, 1>& Vmr_flag,
                   const std::string& Gas_name,
                   boost::shared_ptr<StateMapping> in_map = boost::make_shared<StateMappingLinear>());
  AbsorberVmrLevel(const boost::shared_ptr<Pressure>& Press,
                   const blitz::Array<double, 1>& Vmr,
                   const bool Vmr_flag,
                   const std::string& Gas_name,
                   boost::shared_ptr<StateMapping> in_map = boost::make_shared<StateMappingLinear>());
  virtual boost::shared_ptr<AbsorberVmr> clone() const;
  %python_attribute(sub_state_identifier, std::string);
  virtual std::string state_vector_name_i(int i) const;
  %python_attribute(pressure_profile, blitz::Array<double, 1>);
  %python_attribute(vmr_profile, blitz::Array<double, 1>);
  %python_attribute(vmr_covariance, blitz::Array<double, 2>);
  %python_attribute(vmr_uncertainty, blitz::Array<double, 1>);
  %pickle_serialization();
protected:
  virtual void calc_vmr() const;
};
}


// -*- mode: c++; -*-
// (Not really c++, but closest emacs mode)

%include "fp_common.i"
%{
#include "absorber_vmr_scaled.h"
%}
%base_import(absorber_vmr_imp_base)

%fp_shared_ptr(FullPhysics::AbsorberVmrScaled)
namespace FullPhysics {
class AbsorberVmrScaled : public AbsorberVmrImpBase {
public:
  AbsorberVmrScaled(const boost::shared_ptr<Pressure>& Press,
                    double Scale,                         
                    const std::string& Gas_name);
  virtual boost::shared_ptr<AbsorberVmr> clone() const = 0;
  %python_attribute(sub_state_identifier, std::string);
  virtual std::string state_vector_name_i(int i) const;
  %python_attribute(scale_factor, double)
  %python_attribute(scale_uncertainty, double)
  %pickle_serialization();
protected:
  virtual void calc_vmr() const;
};
}


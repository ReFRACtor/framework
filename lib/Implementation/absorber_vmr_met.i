// -*- mode: c++; -*-
// (Not really c++, but closest emacs mode)
%include "fp_common.i"
%{
#include "absorber_vmr_met.h"
%}
%base_import(absorber_vmr_scaled)
%import "meteorology.i"
%import "pressure.i"
%fp_shared_ptr(FullPhysics::AbsorberVmrMet)
namespace FullPhysics {

%feature("notabstract") AbsorberVmrMet;

class AbsorberVmrMet : public AbsorberVmrScaled {
public:
  AbsorberVmrMet(const boost::shared_ptr<Meteorology>& Met_file,
		   const boost::shared_ptr<Pressure>& Press,
		   double Scale,                         
		   const std::string& Gas_name);
  virtual boost::shared_ptr<AbsorberVmr> clone() const;
  %python_attribute(scale_factor, double)
  %python_attribute(scale_uncertainty, double)
};
}


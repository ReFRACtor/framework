// -*- mode: c++; -*-
// (Not really c++, but closest emacs mode)

%include "fp_common.i"

%{
#include "optical_properties_wrt_rt.h"
#include "sub_state_vector_array.h"
#include "altitude.h"
%}

%base_import(optical_properties_init_base)

%fp_shared_ptr(FullPhysics::OpticalPropertiesWrtRt)

%template (vector_optical_properties_wrt_rt) std::vector<boost::shared_ptr<FullPhysics::OpticalPropertiesWrtRt> >;

namespace FullPhysics {
class OpticalPropertiesWrtRt : public OpticalPropertiesInitBase {
public:
  OpticalPropertiesWrtRt() : OpticalPropertiesInitBase() {};
  %pickle_serialization();
protected:
  virtual void initialize_with_jacobians(const ArrayAd<double, 1>& rayleigh_od, 
					 const ArrayAd<double, 2>& gas_od,
					 const ArrayAd<double, 2>& aerosol_ext_od,
					 const ArrayAd<double, 2>& aerosol_sca_od,
					 const boost::shared_ptr<AerosolPhaseFunctionHelper>& aer_pf_helper);
};

}


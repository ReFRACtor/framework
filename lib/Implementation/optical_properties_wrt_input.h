#ifndef OPTICAL_PROP_WRT_INPUT_H
#define OPTICAL_PROP_WRT_INPUT_H

#include "optical_properties_init_base.h"

namespace FullPhysics {

/****************************************************************//**
  Represents the optical properties for a single spectral point where 
  jacobians are with respect to the input jacobians. Therefore the
  value of the intermediate_jacobian is an identity matrix.
 *******************************************************************/

class OpticalPropertiesWrtInput : public virtual OpticalPropertiesInitBase {
public:
  OpticalPropertiesWrtInput() : OpticalPropertiesInitBase() {};
protected:
  virtual void initialize_with_jacobians(const ArrayAd<double, 1>& rayleigh_od, 
					 const ArrayAd<double, 2>& gas_od,
					 const ArrayAd<double, 2>& aerosol_ext_od,
					 const ArrayAd<double, 2>& aerosol_sca_od,
					 const boost::shared_ptr<AerosolPhaseFunctionHelper>& aer_pf_helper);
private:
  friend class boost::serialization::access;
  template<class Archive>
  void serialize(Archive & ar, const unsigned int version);
};

}

FP_EXPORT_KEY(OpticalPropertiesWrtInput);
#endif

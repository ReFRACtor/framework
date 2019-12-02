#ifndef PCA_OPTICAL_PROP_WRT_INPUT_H
#define PCA_OPTICAL_PROP_WRT_INPUT_H

#include "optical_properties_imp_base.h"

namespace FullPhysics {

/****************************************************************//**
  Represents the optical properties for a single spectral point where 
  jacobians are with respect to the input jacobians. Therefore the
  value of the intermediate_jacobian is an identity matrix.
 *******************************************************************/

class OpticalPropertiesWrtInput : public virtual OpticalPropertiesImpBase {
public:

    OpticalPropertiesWrtInput() : OpticalPropertiesImpBase() {};

protected:

    virtual void initialize_with_jacobians(const ArrayAd<double, 1>& rayleigh_od, 
                                           const ArrayAd<double, 2>& gas_od,
                                           const ArrayAd<double, 2>& aerosol_ext_od,
                                           const ArrayAd<double, 2>& aerosol_sca_od,
                                           const std::vector<ArrayAd<double, 3> >& aerosol_pf_moments);

};

}

#endif

#ifndef OPTICAL_PROP_WRT_RT_H
#define OPTICAL_PROP_WRT_RT_H

#include "optical_properties_imp_base.h"

namespace FullPhysics {

/****************************************************************//**
  Represents the optical properties for a single spectral point where 
  jacobians are with respect to the radiative transfer parameters,
  specifically:
  * gas optical depth per layer (sum of all particles)
  * rayleigh optical depth
  * aerosol optical depth per particle

  The different jacobian basis helps improve the speed of radiative
  transfer computations which are generally linear in speed with the
  number of weigthing functions used. After radiative transfer computation
  the jacobians with respect to the input parameters can be obtained by
  multiplying by the intermediate_jacobian value.
  *******************************************************************/

class OpticalPropertiesWrtRt : public virtual OpticalPropertiesImpBase {
public:

    OpticalPropertiesWrtRt() : OpticalPropertiesImpBase() {};

protected:

    virtual void initialize_with_jacobians(const ArrayAd<double, 1>& rayleigh_od, 
                                           const ArrayAd<double, 2>& gas_od,
                                           const ArrayAd<double, 2>& aerosol_ext_od,
                                           const ArrayAd<double, 2>& aerosol_sca_od,
                                           const boost::shared_ptr<AerosolPhaseFunctionHelper>& aer_pf_helper);

    virtual const std::vector<ArrayAd<double, 3> > aerosol_phase_function_moments_per_particle(int num_moments = -1, int num_scattering = -1) const;

};

}

#endif

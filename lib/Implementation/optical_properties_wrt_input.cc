#include "optical_properties_wrt_input.h"
#include "fp_serialize_support.h"

using namespace blitz;
using namespace FullPhysics;

#ifdef FP_HAVE_BOOST_SERIALIZATION
template<class Archive>
void OpticalPropertiesWrtInput::serialize(Archive & ar,
			const unsigned int UNUSED(version))
{
  ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(OpticalPropertiesInitBase);
}

FP_IMPLEMENT(OpticalPropertiesWrtInput);
#endif

//-----------------------------------------------------------------------
/// Initialize internal variables along with jacobians.
///
/// Jacobians are with respect to the input
/// variables and therefore the intermdiate jacobian is just an identity 
/// matrix.
//-----------------------------------------------------------------------

void OpticalPropertiesWrtInput::initialize_with_jacobians(const ArrayAd<double, 1>& rayleigh_od, 
                                                  const ArrayAd<double, 2>& gas_od,
                                                  const ArrayAd<double, 2>& aerosol_ext_od,
                                                  const ArrayAd<double, 2>& aerosol_sca_od,
                                                  const boost::shared_ptr<AerosolPhaseFunctionHelper>& aer_pf_helper)
{

    // Take references instead of copying for speed
    rayleigh_optical_depth_.reference(rayleigh_od);
    gas_optical_depth_per_particle_.reference(gas_od);
    aerosol_extinction_optical_depth_per_particle_.reference(aerosol_ext_od);
    aerosol_scattering_optical_depth_per_particle_.reference(aerosol_sca_od);

    aerosol_phase_function_helper_ = aer_pf_helper;

    int num_layers = rayleigh_optical_depth_.rows();
    int num_jac = rayleigh_optical_depth_.number_variable();

    // Interermediate jacobian here is the identity matrix
    intermediate_jacobian_.resize(num_layers, num_jac, num_jac);
    intermediate_jacobian_ = 0.0;
    for(int lay_idx = 0; lay_idx < rayleigh_optical_depth_.rows(); lay_idx++) {
        for(int jac_idx = 0; jac_idx < rayleigh_optical_depth_.number_variable(); jac_idx++) {
            intermediate_jacobian_(lay_idx, jac_idx, jac_idx) = 1.0;
        }
    }

} 

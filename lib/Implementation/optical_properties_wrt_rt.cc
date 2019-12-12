#include "optical_properties_wrt_rt.h"

using namespace blitz;
using namespace FullPhysics;

// Index of parameters in jacobians for OpticalPropertiesWrtRt
// Order is due to heritage
const int gas_jac_index = 0;
const int rayleigh_jac_index = 1;
const int aerosol_0_jac_index = 2;

//-----------------------------------------------------------------------
/// Initialize internal variables along with jacobians.
///
/// Jacobians are with respect to the following parameters in this order:
///  * gas optical depth per layer (sum of all particles)
///  * rayleigh optical depth
///  * aerosol optical depth per particle
////-----------------------------------------------------------------------

void OpticalPropertiesWrtRt::initialize_with_jacobians(const ArrayAd<double, 1>& rayleigh_od, 
                                                       const ArrayAd<double, 2>& gas_od,
                                                       const ArrayAd<double, 2>& aerosol_ext_od,
                                                       const ArrayAd<double, 2>& aerosol_sca_od,
                                                       const boost::shared_ptr<AerosolPhaseFunctionHelper>& aer_pf_helper)
{

    Range ra = Range::all();
    firstIndex i1; secondIndex i2;

    // Jacobians for gas optical depth, rayleigh optical depth and aerosols
    int num_layers = rayleigh_od.rows();
    int num_gas = gas_od.cols();
    int num_aer = aerosol_ext_od.cols();

    int num_inp_jac = rayleigh_od.number_variable();
    int num_interm_jac = 2 + num_aer;

    // Create mapping of intermediate jacobians to input jacobians for each layer
    intermediate_jacobian_.resize(num_layers, num_interm_jac, num_inp_jac);

    // At gas optical depth intermediate index the jacobian is the sum of all gas particles, therefore jacobian
    // for each particle is a fraction of the total
    gas_optical_depth_per_particle_.value().reference(gas_od.value());

    Array<double, 3> gas_jac(num_layers, num_gas, num_interm_jac);
    gas_jac = 0.0;
    for(int gas_idx = 0; gas_idx < num_gas; gas_idx++) {
        gas_jac(ra, gas_idx, gas_jac_index) = 1.0/num_gas;
    }
    gas_optical_depth_per_particle_.jacobian().reference(gas_jac);

    for(int lay_idx = 0; lay_idx < gas_od.rows(); lay_idx++) {
        intermediate_jacobian_(lay_idx, gas_jac_index, ra) =  sum(gas_od.jacobian()(lay_idx, ra, ra)(i2, i1), i2);
    }

    // At rayleigh intermediate index the jacobian is the rayleigh jacobian itself
    rayleigh_optical_depth_.value().reference(rayleigh_od.value());

    Array<double, 2> ray_jac(num_layers, num_interm_jac);
    ray_jac = 0.0;
    ray_jac(ra, rayleigh_jac_index) = 1.0;
    rayleigh_optical_depth_.jacobian().reference(ray_jac);

    intermediate_jacobian_(ra, rayleigh_jac_index, ra) = rayleigh_od.jacobian();

    // At each aerosol intermediate index the jacobian is the aerosol extinction optical depth jacobian
    // the scattering jacobian is not included to reduce the number of jacobians even if it reduces the
    // accuracy of the aerosol jacobians
    aerosol_extinction_optical_depth_per_particle_.value().reference(aerosol_ext_od.value());
    aerosol_scattering_optical_depth_per_particle_.value().reference(aerosol_sca_od.value());

    Array<double, 3> aer_ext_jac(num_layers, num_aer, num_interm_jac);
    Array<double, 3> aer_sca_jac(num_layers, num_aer, num_interm_jac);

    aer_ext_jac = 0.0;
    aer_sca_jac = 0.0;

    for(int aer_idx = 0; aer_idx < num_aer; aer_idx++) {
        aer_ext_jac(ra, aer_idx, aerosol_0_jac_index+aer_idx) = 1.0;

        // This is a consequence of the relationship between tau_ext and tau_sca:
        // tau_sca = tau_ext / ssa_aer
        // Where ssa_aer = k_sca/k_ext
        // k here stand for aerosol coefficient
        // Since tau_sca = k_sca and tau_ext = k_sca (with a wavelength interolation term not shown)
        // We can through the chain rule show that
        // tau_sca' = tau_sca / tau_ext
        // And hence we form a relationship between tau_ext and tau_sca despite tau_sca not being included in 
        // the intermediate jacobian.
        if(!aerosol_ext_od.is_constant() && !aerosol_sca_od.is_constant()) {
            aer_sca_jac(ra, aer_idx, aerosol_0_jac_index+aer_idx) = aerosol_sca_od.value()(ra, aer_idx) / aerosol_ext_od.value()(ra, aer_idx);
        } else {
            aer_sca_jac(ra, aer_idx, aerosol_0_jac_index+aer_idx) = 0.0;
        }

        intermediate_jacobian_(ra, aerosol_0_jac_index+aer_idx, ra) = aerosol_ext_od.jacobian()(ra, aer_idx, ra);
    }

    aerosol_extinction_optical_depth_per_particle_.jacobian().reference(aer_ext_jac);
    aerosol_scattering_optical_depth_per_particle_.jacobian().reference(aer_sca_jac);

    aerosol_phase_function_helper_ = aer_pf_helper;
}

//-----------------------------------------------------------------------
/// Override the default behavior of the parent class so we can ensure
/// that the number of jacobians for this value is indeed constant. For
/// now we are not able to handle that case. But should the need arise
/// this method serves as a placeholder where the jacobian manipulation
/// would need to occur.
//-----------------------------------------------------------------------

const std::vector<ArrayAd<double, 3> > OpticalPropertiesWrtRt::aerosol_phase_function_moments_per_particle(int num_moments, int num_scattering) const
{
    std::vector<ArrayAd<double, 3> > aerosol_pf_moments = OpticalPropertiesImpBase::aerosol_phase_function_moments_per_particle(num_moments, num_scattering);

    for(int part_idx = 0; part_idx < aerosol_pf_moments.size(); part_idx++) {
        if (!aerosol_pf_moments[part_idx].is_constant()) {
            throw Exception("Phase function moments with a jacobian component are not handled at this time");
        }
    }

    return aerosol_pf_moments;
}

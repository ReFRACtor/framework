#include "optical_properties_wrt_rt.h"
#include "fp_serialize_support.h"

using namespace blitz;
using namespace FullPhysics;

#ifdef FP_HAVE_BOOST_SERIALIZATION
template<class Archive>
void OpticalPropertiesWrtRt::serialize(Archive & ar,
                                       const unsigned int UNUSED(version))
{
    ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(OpticalPropertiesInitBase);
}

FP_IMPLEMENT(OpticalPropertiesWrtRt);
#endif

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
                                                       const boost::shared_ptr<AerosolPhaseFunctionHelper>& aer_pf_helper,
                                                       const int num_jacobians)
{

    Range ra = Range::all();
    firstIndex i1;
    secondIndex i2;

    // Jacobians for gas optical depth, rayleigh optical depth and aerosols
    int num_layers = rayleigh_od.rows();
    int num_gas = gas_od.cols();
    int num_aer = aerosol_ext_od.cols();

    // If left unspecified try and get number of jacobians from input objects, but
    // will only work if there is pressure based quantity in the state vector
    int num_inp_jac = num_jacobians;
    if (num_inp_jac <= 0) {
        num_inp_jac = rayleigh_od.number_variable();

        if(num_gas > 0) {
            num_inp_jac = std::max(num_inp_jac, gas_od.number_variable());
        }

        if(num_aer > 0) {
            num_inp_jac = std::max(num_inp_jac, aerosol_ext_od.number_variable());
        }
    }

    if(num_inp_jac > 0 &&
       ((!rayleigh_od.is_constant() && rayleigh_od.number_variable() != num_inp_jac) ||
        (num_gas > 0 && !gas_od.is_constant() && gas_od.number_variable() != num_inp_jac) ||
        (num_aer > 0 && !aerosol_ext_od.is_constant() && aerosol_ext_od.number_variable() != num_inp_jac))) {
        Exception err;
        err << "rayleigh_od (" << rayleigh_od.number_variable() << "), "
            << "gas_od(" << gas_od.number_variable() << ") and "
            << "aerosol_ext_od (" << aerosol_ext_od.number_variable() << ") "
            <<" all need to have the same size jacobian (" << num_inp_jac << ") "
            << "if they are not constant";
        throw err;
    }

    int num_interm_jac = 2 + num_aer;

    // Create mapping of intermediate jacobians to input jacobians for
    // each layer
    intermediate_jacobian_.resize(num_layers, num_interm_jac, num_inp_jac);

    // Gas optical depth per particle ArrayAd will be set to be a constant
    // and hence have no jacobian portion, this is due to the fact that
    // the intermediate jacobian is defined in terms of the total gas optical depth
    // There is no way to compute the actual jacobian for the per particle
    // values that can be transformed by the intermediate jacobian matrix
    //
    // The resize below is to enforce that is_const in the ArrayAd is set to True
    gas_optical_depth_per_particle_.resize(gas_od.value().shape(), 0);
    gas_optical_depth_per_particle_.value().reference(gas_od.value());

    // Go ahead and compute the total gas optical depth values per layer since
    // we will be setting the jacobian values explicitly since they do not flow
    // from the per particle jacobians since they are constant
    gas_optical_depth_per_layer_.resize(gas_od.value().rows(), num_interm_jac);

    if(gas_optical_depth_per_layer_.rows() > 0) {
        gas_optical_depth_per_layer_.value() = sum(gas_od.value()(i1, i2), i2);

        gas_optical_depth_per_layer_.jacobian() = 0.0;
        gas_optical_depth_per_layer_.jacobian()(ra, gas_jac_index) = 1.0;
    }

    // Intermediate value for the total gas optical depth is the sum
    // of the per gas jacobians,
    for(int lay_idx = 0; lay_idx < gas_od.rows(); lay_idx++) {
        if(!gas_od.is_constant()) {
            intermediate_jacobian_(lay_idx, gas_jac_index, ra) =
                sum(gas_od.jacobian()(lay_idx, ra, ra)(i2, i1), i2);
        } else {
            if(num_inp_jac != 0) {
                intermediate_jacobian_(lay_idx, gas_jac_index, ra) = 0;
            }
        }
    }

    // At rayleigh intermediate index the jacobian is the rayleigh jacobian itself
    rayleigh_optical_depth_.value().reference(rayleigh_od.value());

    Array<double, 2> ray_jac(num_layers, num_interm_jac);
    ray_jac = 0.0;
    ray_jac(ra, rayleigh_jac_index) = 1.0;
    rayleigh_optical_depth_.jacobian().reference(ray_jac);

    if(!rayleigh_od.is_constant()) {
        intermediate_jacobian_(ra, rayleigh_jac_index, ra) = rayleigh_od.jacobian();
    } else {
        if(num_inp_jac != 0) {
            intermediate_jacobian_(ra, rayleigh_jac_index, ra) = 0;
        }
    }

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
        aer_ext_jac(ra, aer_idx, aerosol_0_jac_index + aer_idx) = 1.0;

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
            aer_sca_jac(ra, aer_idx, aerosol_0_jac_index + aer_idx) =
                where(aerosol_ext_od.value()(ra, aer_idx) != 0,
                      aerosol_sca_od.value()(ra, aer_idx) / aerosol_ext_od.value()(ra, aer_idx), 0.0);
        } else {
            aer_sca_jac(ra, aer_idx, aerosol_0_jac_index + aer_idx) = 0.0;
        }

        if(!aerosol_ext_od.is_constant()) {
            intermediate_jacobian_(ra, aerosol_0_jac_index + aer_idx, ra) = aerosol_ext_od.jacobian()(ra, aer_idx, ra);
        } else {
            if(num_inp_jac != 0) {
                intermediate_jacobian_(ra, aerosol_0_jac_index + aer_idx, ra) = 0;
            }
        }
    }

    aerosol_extinction_optical_depth_per_particle_.jacobian().reference(aer_ext_jac);
    aerosol_scattering_optical_depth_per_particle_.jacobian().reference(aer_sca_jac);

    aerosol_phase_function_helper_ = aer_pf_helper;
}

//-----------------------------------------------------------------------
/// Override the default behavior. This value is computed in the
/// initialization routine because it is one of the items by which
/// the intermediate jacobian is defined. Since it has already been
/// computed, there is no need to attempt to compute it like in the
/// base class and here we simply return the value created in the
/// init routine.
//-----------------------------------------------------------------------

ArrayAd<double, 1> OpticalPropertiesWrtRt::gas_optical_depth_per_layer() const
{
    return gas_optical_depth_per_layer_;
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

    for(int part_idx = 0; part_idx < (int) aerosol_pf_moments.size(); part_idx++) {
        if (!aerosol_pf_moments[part_idx].is_constant()) {
            throw Exception("Phase function moments with a jacobian component are not handled at this time");
        }
    }

    return aerosol_pf_moments;
}

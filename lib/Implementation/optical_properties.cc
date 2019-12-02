#include "optical_properties.h"

#include "rayleigh_greek_moment.h"

using namespace blitz;
using namespace FullPhysics;

//-----------------------------------------------------------------------
/// Initializes an OpticalProperties object using precomputed values:
/// In the parameters \f$l\f$ stands for layer index, \f$p\f$ stands for particle index.
/// \param rayleigh_od Rayleigh optical depth for each layer (\f$\tau_{ray,l}\f$)
/// \param gas_od Gas Absorber optical depth for each layer and gas particle type (\f$\tau_{gas,lp}\f$)
/// \param aerosol_ext_od Aerosol extinction optical depth for each layer and aerosol particle type (\f$\tau_{aer\_ext,lp}\f$)
/// \param aerosol_sca_od Aerosol scattering optical depth for each layer and aerosol particle type (\f$\tau_{aer\_sca,lp}\f$)
//-----------------------------------------------------------------------

void OpticalPropertiesImpBase::initialize(const ArrayAd<double, 1>& rayleigh_od, 
                                   const ArrayAd<double, 2>& gas_od,
                                   const ArrayAd<double, 2>& aerosol_ext_od,
                                   const ArrayAd<double, 2>& aerosol_sca_od,
                                   const std::vector<ArrayAd<double, 3> >& aerosol_pf_moments)
{
    initialize_with_jacobians(rayleigh_od, gas_od, aerosol_ext_od, aerosol_sca_od, aerosol_pf_moments);

    assert_sizes();
    initialized = true;
}

//-----------------------------------------------------------------------
/// Initializes an OpticalProperties object from atmospheric objects
/// at a specific spectral point for a given instrument channel.
//-----------------------------------------------------------------------

void OpticalPropertiesImpBase::initialize(const DoubleWithUnit spectral_point,
                                   const int channel_index,
                                   const boost::shared_ptr<Absorber>& absorber,
                                   const boost::shared_ptr<Rayleigh>& rayleigh,
                                   const boost::shared_ptr<Aerosol>& aerosol,
                                   int num_pf_mom, int num_scattering)
{

    Range ra = Range::all();
    double wn = spectral_point.convert_wave(units::inv_cm).value;

    ArrayAd<double, 1> rayleigh_od(rayleigh->optical_depth_each_layer(wn, channel_index));
    ArrayAd<double, 2> gas_od(absorber->optical_depth_each_layer(wn, channel_index));
    ArrayAd<double, 2> aerosol_ext_od(aerosol->extinction_optical_depth_each_layer(wn));
    ArrayAd<double, 2> aerosol_sca_od(aerosol->scattering_optical_depth_each_layer(wn));

    // Copy out aerosol phase function moments
    std::vector<ArrayAd<double, 3> > aerosol_pf_moments;
    for(int p_idx = 0; p_idx < aerosol->number_particle(); ++p_idx) {
        aerosol_pf_moments.push_back(aerosol->pf_mom(wn, p_idx, num_pf_mom, num_scattering));
    }

    initialize_with_jacobians(rayleigh_od, gas_od, aerosol_ext_od, aerosol_sca_od, aerosol_pf_moments);

    assert_sizes();
    initialized = true;

}

//-----------------------------------------------------------------------
/// Protected method that throws an exception if the class has not yet
/// been initialized. Due to the desire to allow extension by extending
/// the initialize_with_jacobians method, having a two phase initializtion
/// is necessary. This method checks that the second phase has been 
/// run.
//-----------------------------------------------------------------------

void OpticalPropertiesImpBase::assert_init() const
{
    if (!initialized) {
        throw Exception("This optical properties instance has not yet been initalized. Call an initialize method.");
    }
}

//-----------------------------------------------------------------------
/// Protected method that throws an exception if the size of input
/// variables are not consistent with each other.
//-----------------------------------------------------------------------

void OpticalPropertiesImpBase::assert_sizes() const
{
    int num_layers = rayleigh_optical_depth_.rows();

    if(gas_optical_depth_per_particle_.rows() != num_layers) {
        throw Exception("Gas optical depth value has inconsistent number of layers");
    }

    if (aerosol_extinction_optical_depth_per_particle_.rows() != num_layers) {
        throw Exception("Aerosol extinction optical depth has inconsistent number of layers");
    }
    if (aerosol_scattering_optical_depth_per_particle_.rows() != num_layers) {
        throw Exception("Aerosol scattering optical depth has inconsistent number of layers");
    }

    if (aerosol_extinction_optical_depth_per_particle_.cols() != aerosol_scattering_optical_depth_per_particle_.cols()) {
        throw Exception("Aerosol extinction and scattering optical depth has inconsistent number of particles");
    }

    if (aerosol_extinction_optical_depth_per_particle_.cols() != aerosol_phase_function_moments_per_particle_.size()) {
        throw Exception("Aerosol extinction optical depth and aerosol phase function moments have inconsistent number of particles");
    }
}

//-----------------------------------------------------------------------
/// Gas Absorber optical depth summed over all particles:
/// \f$ \tau_{gas,l} = \sum_{p=0}^{P} \tau_{gas,lp} \f$
//-----------------------------------------------------------------------

ArrayAd<double, 1> OpticalPropertiesImpBase::gas_optical_depth_per_layer() const
{
    assert_init();

    // gas component
    ArrayAd<double, 2> gas_od_per_p(gas_optical_depth_per_particle());
    if(gas_optical_depth_per_layer_.rows() == 0 and gas_od_per_p.cols() > 0) {
        firstIndex i1; secondIndex i2;
        Range ra = Range::all();

        gas_optical_depth_per_layer_.resize(gas_od_per_p.rows(), gas_od_per_p.number_variable());

        gas_optical_depth_per_layer_.value() = sum(gas_od_per_p.value(), i2);
        if (!gas_od_per_p.is_constant()) {
            for(int lay_idx = 0; lay_idx < gas_od_per_p.value().rows(); lay_idx++) {
                gas_optical_depth_per_layer_.jacobian()(lay_idx, ra) = sum(gas_od_per_p.jacobian()(lay_idx, ra, ra)(i2, i1), i2);
            }
        }
    }
    return gas_optical_depth_per_layer_;
}

//-----------------------------------------------------------------------
/// Aerosol extinction optical depth summed over all particles
/// \f$ \tau_{aer\_ext,l} = \sum_{p=0}^{P} \tau_{aer\_ext,lp} \f$
//-----------------------------------------------------------------------

ArrayAd<double, 1> OpticalPropertiesImpBase::aerosol_extinction_optical_depth_per_layer() const
{
    assert_init(); 

    // aerosol component
    ArrayAd<double, 2> aer_od_per_p(aerosol_extinction_optical_depth_per_particle());
    if (aerosol_extinction_optical_depth_per_layer_.rows() == 0 and aer_od_per_p.cols() > 0) {
        firstIndex i1; secondIndex i2;
        Range ra = Range::all();

        aerosol_extinction_optical_depth_per_layer_.resize(aer_od_per_p.rows(), aer_od_per_p.number_variable());

        aerosol_extinction_optical_depth_per_layer_.value() = sum(aer_od_per_p.value(), i2);
        if (!aer_od_per_p.is_constant()) {
            for(int lay_idx = 0; lay_idx < aer_od_per_p.value().rows(); lay_idx++) {
                aerosol_extinction_optical_depth_per_layer_.jacobian()(lay_idx, ra) = sum(aer_od_per_p.jacobian()(lay_idx, ra, ra)(i2, i1), i2);
            }
        }
    }

    return aerosol_extinction_optical_depth_per_layer_;
}


//-----------------------------------------------------------------------
/// Aerosol scattering optical depth summed over all particles
/// \f$ \tau_{aer\_sca,l} = \sum_{p=0}^{P} \tau_{aer\_sca,lp} \f$
//-----------------------------------------------------------------------

ArrayAd<double, 1> OpticalPropertiesImpBase::aerosol_scattering_optical_depth_per_layer() const
{
    assert_init(); 

    ArrayAd<double, 2> aer_ssa_per_p(aerosol_scattering_optical_depth_per_particle());
    if (aerosol_scattering_optical_depth_per_layer_.rows() == 0 and aer_ssa_per_p.cols() > 0) {
        firstIndex i1; secondIndex i2;
        Range ra = Range::all();

        aerosol_scattering_optical_depth_per_layer_.resize(aer_ssa_per_p.rows(), aer_ssa_per_p.number_variable());

        aerosol_scattering_optical_depth_per_layer_.value() = sum(aer_ssa_per_p.value(), i2);
        if (!aer_ssa_per_p.is_constant()) {
            for(int lay_idx = 0; lay_idx < aer_ssa_per_p.value().rows(); lay_idx++) {
                aerosol_scattering_optical_depth_per_layer_.jacobian()(lay_idx, ra) = sum(aer_ssa_per_p.jacobian()(lay_idx, ra, ra)(i2, i1), i2);
            }
        }
    }

    return aerosol_scattering_optical_depth_per_layer_;
}

//-----------------------------------------------------------------------
/// Total optical depth consisting of rayleigh + gas + aerosol
/// Gas and aerosol contributions are allowed to possibly be zero.
/// Computed as:
/// \f$ \tau_{tot,l} = \tau_{ray,l} + \tau_{gas,l}  + \tau_{aer\_ext,l} \f$
//-----------------------------------------------------------------------

ArrayAd<double, 1> OpticalPropertiesImpBase::total_optical_depth() const
{
    assert_init(); 

    if (total_optical_depth_.rows() == 0) {

        total_optical_depth_.resize(rayleigh_optical_depth_.rows(), rayleigh_optical_depth_.number_variable());

        // rayleigh component
        // Use rayleigh to initalize the total_optical_depth values
        ArrayAd<double, 1> rayleigh_od(rayleigh_optical_depth());
        total_optical_depth_.value() = rayleigh_od.value();
        if(!rayleigh_od.is_constant()) {
            total_optical_depth_.jacobian() = rayleigh_od.jacobian();
        } else {
            total_optical_depth_.jacobian() = 0;
        }

        ArrayAd<double, 1> gas_od(gas_optical_depth_per_layer());
        if(gas_od.rows() > 0) {
            total_optical_depth_.value() += gas_od.value();
            if(!gas_od.is_constant()) {
                total_optical_depth_.jacobian() += gas_od.jacobian();
            }
        }
 
        ArrayAd<double, 1> aerosol_od(aerosol_extinction_optical_depth_per_layer());
        if(aerosol_od.rows() > 0) {
            total_optical_depth_.value() += aerosol_od.value();
            if(!aerosol_od.is_constant()) {
                total_optical_depth_.jacobian() += aerosol_od.jacobian();
            }
        }
 
    }
    
    return total_optical_depth_;
}

// Private helper that just computes the sum of the rayleigh and aerosol
// optical depths, the value is used in multiple other calculations and
// computed once for speed purposes
ArrayAd<double, 1> OpticalPropertiesImpBase::scattering_sum() const
{
    assert_init(); 

    if(scattering_sum_.rows() == 0) {
        ArrayAd<double, 1> ray_od(rayleigh_optical_depth());

        scattering_sum_.resize(ray_od.rows(), ray_od.number_variable());

        scattering_sum_.value() = ray_od.value();
        if(!ray_od.is_constant()) {
            scattering_sum_.jacobian() = ray_od.jacobian();
        }

        ArrayAd<double, 1> aer_sca_od(aerosol_scattering_optical_depth_per_layer());
        if(aer_sca_od.rows() > 0) {
            scattering_sum_.value() += aer_sca_od.value();
            if(!aer_sca_od.is_constant()) {
                scattering_sum_.jacobian() += aer_sca_od.jacobian();
            }
        }

    }
    
    return scattering_sum_;
}

//-----------------------------------------------------------------------
/// Total single scattering albedo
/// Aerosol contribution is allowed to possibly be zero.
/// Computed as:
/// \f$ \omega_{tot,l} = \frac{ \tau_{ray,l} + \tau_{aer\_sca,l} }{ \tau_{tot,l} } \f$
//-----------------------------------------------------------------------

ArrayAd<double, 1> OpticalPropertiesImpBase::total_single_scattering_albedo() const
{
    assert_init(); 

    if(total_single_scattering_albedo_.rows() == 0) {
        firstIndex i1; 

        ArrayAd<double, 1> ssa_sum(scattering_sum());
        ArrayAd<double, 1> tot_od(total_optical_depth());

        total_single_scattering_albedo_.resize(ssa_sum.rows(), ssa_sum.number_variable());

        total_single_scattering_albedo_.value() = ssa_sum.value() / tot_od.value();

        if(!total_single_scattering_albedo_.is_constant()) {
            // Jacobian computed using quotient to avoid divide by zero possible in tot_od jacobian
            total_single_scattering_albedo_.jacobian() = 
                ssa_sum.jacobian() / tot_od.value()(i1) - 
                ssa_sum.value()(i1) / (tot_od.value()(i1) * tot_od.value()(i1)) * tot_od.jacobian();
        }

    }

    return total_single_scattering_albedo_;
}

//-----------------------------------------------------------------------
/// Fraction of the scattering attributable to rayleigh at each layer
//-----------------------------------------------------------------------

ArrayAd<double, 1> OpticalPropertiesImpBase::rayleigh_fraction() const
{
    assert_init(); 

    if(rayleigh_fraction_.rows() == 0) {
        firstIndex i1; 

        ArrayAd<double, 1> ray_od(rayleigh_optical_depth());
        ArrayAd<double, 1> ssa_sum(scattering_sum());

        rayleigh_fraction_.resize(ray_od.rows(), ray_od.number_variable());
        rayleigh_fraction_.value() = ray_od.value() / ssa_sum.value();

        if(!rayleigh_fraction_.is_constant()) {
            // Jacobian computed using quotient to avoid divide by zero possible in tot_od jacobian
            rayleigh_fraction_.jacobian() = 
                ray_od.jacobian() / ssa_sum.value()(i1) - 
                ray_od.value()(i1) / (ssa_sum.value()(i1) * ssa_sum.value()(i1)) * ssa_sum.jacobian();
        }
    }

    return rayleigh_fraction_;
}

//-----------------------------------------------------------------------
/// Fraction of the scattering attributable to each aerosol particle at each layer
//-----------------------------------------------------------------------

ArrayAd<double, 2> OpticalPropertiesImpBase::aerosol_fraction() const
{
    assert_init(); 

    if(aerosol_fraction_.rows() == 0) {
        firstIndex i1; secondIndex i2; thirdIndex i3;

        ArrayAd<double, 2> aer_sca(aerosol_scattering_optical_depth_per_particle());
        ArrayAd<double, 1> ssa_sum(scattering_sum());

        aerosol_fraction_.resize(aer_sca.rows(), aer_sca.cols(), aer_sca.number_variable());
        aerosol_fraction_.value() = aer_sca.value() / ssa_sum.value()(i1);

        if(!aerosol_fraction_.is_constant()) {
            // Jacobian computed using quotient to avoid divide by zero possible in tot_od jacobian
            aerosol_fraction_.jacobian() = 
                aer_sca.jacobian() / ssa_sum.value()(i1) - 
                aer_sca.value()(i1,i2) / (ssa_sum.value()(i1) * ssa_sum.value()(i1)) * ssa_sum.jacobian()(i1, i3);
        }
    }

    return aerosol_fraction_;
}


//-----------------------------------------------------------------------
/// Returns the portion of the total phase function attributable to
/// aerosols. This is the summation over particles multiplied by the
/// fraction of aerosol.
///
/// Output dimensions: num_moments x num_layers x num_scattering
//-----------------------------------------------------------------------

ArrayAd<double, 3> OpticalPropertiesImpBase::aerosol_phase_function_moments_portion() const
{
    firstIndex i1; secondIndex i2; thirdIndex i3; fourthIndex i4;
    Range ra = Range::all();

    if(aerosol_phase_function_moments_portion_.rows() == 0 and aerosol_phase_function_moments_per_particle_.size() > 0) {
        // Gather the maximum number of moments and scattering dimensions

        int max_num_moms = 0; int max_num_scatt = 0; int num_var = -1;
        for(int part_idx = 0; part_idx < aerosol_phase_function_moments_per_particle_.size(); part_idx++) {
            max_num_moms = std::max(max_num_moms, aerosol_phase_function_moments_per_particle_[part_idx].rows());
            max_num_scatt = std::max(max_num_scatt, aerosol_phase_function_moments_per_particle_[part_idx].depth());

            if(num_var > 0 && num_var != aerosol_phase_function_moments_per_particle_[part_idx].number_variable()) {
                throw Exception("Phase function moments per particle have inconsistent number of jacobians");
            }
            num_var = aerosol_phase_function_moments_per_particle_[part_idx].number_variable();
        }

        ArrayAd<double, 2> frac_aer(aerosol_fraction());

        if (num_var == 0 && !frac_aer.is_constant()) {
            num_var = frac_aer.number_variable();
        }

        aerosol_phase_function_moments_portion_.resize(max_num_moms, number_layers(), max_num_scatt, num_var);
        aerosol_phase_function_moments_portion_.value() = 0;

        if (!aerosol_phase_function_moments_portion_.is_constant()) {
            aerosol_phase_function_moments_portion_.jacobian() = 0;
        }

        for(int part_idx = 0; part_idx < aerosol_phase_function_moments_per_particle_.size(); part_idx++) {
            ArrayAd<double, 3> pf_mom_part(aerosol_phase_function_moments_per_particle_[part_idx]);

            Range r_mom = Range(0, pf_mom_part.rows() - 1);
            Range r_scatt = Range(0, pf_mom_part.depth() - 1);

            Array<double, 3> pf_per_lay_val(aerosol_phase_function_moments_portion_.value()(r_mom, ra, r_scatt));
            Array<double, 4> pf_per_lay_jac(aerosol_phase_function_moments_portion_.jacobian()(r_mom, ra, r_scatt, ra));

            Array<double, 1> fa_val(frac_aer.value()(ra, part_idx));
            Array<double, 2> fa_jac(frac_aer.jacobian()(ra, part_idx, ra));

            pf_per_lay_val += fa_val(i2) * pf_mom_part.value()(i1, i2, i3);

            if(pf_mom_part.is_constant()) {
                if(!frac_aer.is_constant()) {
                    pf_per_lay_jac += fa_jac(i2, i4) * pf_mom_part.value()(i1, i2, i3);
                }
            } else {
                if(frac_aer.is_constant()) {
                    pf_per_lay_jac += fa_val(i2) * pf_mom_part.jacobian()(i1, i2, i3, i4);
                } else {
                    pf_per_lay_jac += fa_jac(i2,i4) * pf_mom_part.value()(i1, i2, i3) + fa_val(i2) * pf_mom_part.jacobian()(i1, i2, i3, i4);
                }
            }

        }

    }

    return aerosol_phase_function_moments_portion_;
}

//-----------------------------------------------------------------------
/// Compute the portion of the total phase function attributable to
/// rayleigh scattering. This returns the Rayleigh Greek Moments matrix
/// multiplied by the fraction of rayleigh.
///
/// Returned dimensions: num_moments x num_layers x num_scattering
//-----------------------------------------------------------------------

ArrayAd<double, 3> OpticalPropertiesImpBase::rayleigh_phase_function_moments_portion() const
{
    if(rayleigh_phase_function_moments_portion_.rows() == 0) {
        firstIndex i1; secondIndex i2; thirdIndex i3; fourthIndex i4;

        Array<double, 2> pf_ray(RayleighGreekMoment::array().copy());
        ArrayAd<double, 1> frac_ray(rayleigh_fraction());

        rayleigh_phase_function_moments_portion_.resize(pf_ray.rows(), number_layers(), pf_ray.cols(), frac_ray.number_variable());

        rayleigh_phase_function_moments_portion_.value() = frac_ray.value()(i2) * pf_ray(i1,i3);
        rayleigh_phase_function_moments_portion_.jacobian() = frac_ray.jacobian()(i2,i4) * pf_ray(i1,i3);
    }

    return rayleigh_phase_function_moments_portion_;
}

//-----------------------------------------------------------------------
/// Compute the total phase function moments which is the summation
/// of the rayleigh and aerosol portions
//-----------------------------------------------------------------------

ArrayAd<double, 3> OpticalPropertiesImpBase::total_phase_function_moments() const
{
    if(total_phase_function_moments_.rows() == 0) {
        ArrayAd<double, 3> pf_ray(rayleigh_phase_function_moments_portion());
        ArrayAd<double, 3> pf_aer(aerosol_phase_function_moments_portion());

        Range r_ray_moms(0, pf_ray.rows() - 1);
        Range ra = Range::all();

        int num_jac = 0;
        if (!pf_aer.is_constant()) {
            num_jac = pf_aer.number_variable();
        } else if(!pf_ray.is_constant()) {
            num_jac = pf_ray.number_variable();
        }

        total_phase_function_moments_.resize(pf_aer.shape(), num_jac);

        total_phase_function_moments_.value() = pf_aer.value();
        total_phase_function_moments_.value()(r_ray_moms, ra, ra) += pf_ray.value();

        total_phase_function_moments_.jacobian() = 0;
        if (!pf_aer.is_constant()) {
            total_phase_function_moments_.jacobian() += pf_aer.jacobian();
        }

        if (!pf_ray.is_constant()) {
            total_phase_function_moments_.jacobian()(r_ray_moms, ra, ra, ra) += pf_ray.jacobian();
        }
    }

    return total_phase_function_moments_;
}

//=======================================================================

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
                                                  const std::vector<ArrayAd<double, 3> >& aerosol_pf_moments)
{

    // Take references instead of copying for speed
    rayleigh_optical_depth_.reference(rayleigh_od);
    gas_optical_depth_per_particle_.reference(gas_od);
    aerosol_extinction_optical_depth_per_particle_.reference(aerosol_ext_od);
    aerosol_scattering_optical_depth_per_particle_.reference(aerosol_sca_od);
    aerosol_phase_function_moments_per_particle_ = aerosol_pf_moments;

    intermediate_jacobian_.resize(rayleigh_optical_depth_.number_variable(), rayleigh_optical_depth_.number_variable());
    intermediate_jacobian_ = 0.0;
    for(int lay_idx = 0; lay_idx < rayleigh_optical_depth_.rows(); lay_idx++) {
        for(int jac_idx = 0; jac_idx < rayleigh_optical_depth_.number_variable(); jac_idx++) {
            intermediate_jacobian_(lay_idx, jac_idx, jac_idx) = 1.0;
        }
    }

} 

//=======================================================================

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
                                                       const std::vector<ArrayAd<double, 3> >& aerosol_pf_moments)
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
        aer_sca_jac(ra, aer_idx, aerosol_0_jac_index+aer_idx) = aerosol_sca_od.value()(ra, aer_idx) / aerosol_ext_od.value()(ra, aer_idx);

        intermediate_jacobian_(ra, aerosol_0_jac_index+aer_idx, ra) = aerosol_ext_od.jacobian()(ra, aer_idx, ra);
    }

    aerosol_extinction_optical_depth_per_particle_.jacobian().reference(aer_ext_jac);
    aerosol_scattering_optical_depth_per_particle_.jacobian().reference(aer_sca_jac);

    aerosol_phase_function_moments_per_particle_ = aerosol_pf_moments;

    for(int part_idx = 0; part_idx < aerosol_pf_moments.size(); part_idx++) {
        if (!aerosol_pf_moments[part_idx].is_constant()) {
            throw Exception("Phase function moments with a jacobian component are not handled at this time");
        }
    }
}

#include "optical_properties.h"

using namespace blitz;
using namespace FullPhysics;

OpticalProperties::OpticalProperties(const ArrayAd<double, 1>& rayleigh_od, 
                                     const ArrayAd<double, 2>& gas_od,
                                     const ArrayAd<double, 2>& aerosol_ext_od,
                                     const ArrayAd<double, 2>& aerosol_sca_od)
: rayleigh_optical_depth_(rayleigh_od),
  gas_optical_depth_per_particle_(gas_od), 
  aerosol_extinction_optical_depth_per_particle_(aerosol_ext_od), 
  aerosol_scattering_optical_depth_per_particle_(aerosol_sca_od)
{
    // The magic is in the initialization
}

OpticalProperties::OpticalProperties(const DoubleWithUnit spectral_point,
                                     const int channel_index,
                                     const boost::shared_ptr<Absorber>& absorber,
                                     const boost::shared_ptr<Rayleigh>& rayleigh,
                                     const boost::shared_ptr<Aerosol>& aerosol)
{
    double wn = spectral_point.convert_wave(units::inv_cm).value;

    rayleigh_optical_depth_ = rayleigh->optical_depth_each_layer(wn, channel_index);
    gas_optical_depth_per_particle_ = absorber->optical_depth_each_layer(wn, channel_index);
    aerosol_extinction_optical_depth_per_particle_ = aerosol->extinction_optical_depth_each_layer(wn);

    aerosol_scattering_optical_depth_per_particle_.resize(aerosol_extinction_optical_depth_per_particle_.shape(), aerosol_extinction_optical_depth_per_particle_.number_variable());
    for(int particle_idx = 0; particle_idx < aerosol_extinction_optical_depth_per_particle_.cols(); ++particle_idx) {
      aerosol_scattering_optical_depth_per_particle_(Range::all(), particle_idx) = aerosol->scattering_optical_depth_each_layer(wn, particle_idx, aerosol_extinction_optical_depth_per_particle_(Range::all(), particle_idx));
    }

}

/// Gas optical depth summed over all particles

ArrayAd<double, 1> OpticalProperties::gas_optical_depth_per_layer() const
{
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

/// Aerosol optical depth summed over all particles

ArrayAd<double, 1> OpticalProperties::aerosol_extinction_optical_depth_per_layer() const
{
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

/// Aerosol single scattering albedo summed over all particles

ArrayAd<double, 1> OpticalProperties::aerosol_scattering_optical_depth_per_layer() const
{
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


/// Total optical depth consisting of rayleigh + gas + aerosol
/// Gas and aerosol contributions may be empty

ArrayAd<double, 1> OpticalProperties::total_optical_depth() const
{
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

ArrayAd<double, 1> OpticalProperties::total_single_scattering_albedo() const
{
    if(total_single_scattering_albedo_.rows() == 0) {

        ArrayAd<double, 1> ray_od(rayleigh_optical_depth());
        ArrayAd<double, 1> tot_od(total_optical_depth());

        total_single_scattering_albedo_.resize(tot_od.rows(), tot_od.number_variable());

        total_single_scattering_albedo_.value() = ray_od.value();
        if(!ray_od.is_constant()) {
            total_single_scattering_albedo_.jacobian() = ray_od.jacobian();
        }

        ArrayAd<double, 1> aer_sca_od(aerosol_scattering_optical_depth_per_layer());
        if(aer_sca_od.rows() > 0) {
            total_single_scattering_albedo_.value() += aer_sca_od.value();
            if(!aer_sca_od.is_constant()) {
                total_single_scattering_albedo_.jacobian() += aer_sca_od.jacobian();
            }
        }

        total_single_scattering_albedo_.value() /= tot_od.value();
        if(!tot_od.is_constant()) {
            total_single_scattering_albedo_.jacobian() /= tot_od.jacobian();
        }

    }

    return total_single_scattering_albedo_;
}

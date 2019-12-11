#ifndef OPTICAL_PROPERTIES_H
#define OPTICAL_PROPERTIES_H

#include "generic_object.h"

#include "array_ad.h"

namespace FullPhysics {

/****************************************************************//**
  Represents the optical properties for a single spectral point.
  The intention is to contain all the information calculated from
  the atmosphere needed for a radiative transfer calculation along
  with intermediate computations useful by approximation methods.
 *******************************************************************/

class OpticalProperties : public virtual GenericObject {
public:

    // Sizes of types of information stored

    virtual int number_layers() const = 0;
    virtual int number_gas_particles() const = 0;
    virtual int number_aerosol_particles() const = 0;
 
    /// Returns Rayleigh optical depth for each layer: \f$ \tau_{ray,l} \f$
    virtual ArrayAd<double, 1> rayleigh_optical_depth() const = 0;

    /// Returns Gas Absorber optical depth for each layer and particle type: \f$ \tau_{gas,lp} \f$
    virtual ArrayAd<double, 2> gas_optical_depth_per_particle() const = 0;

    /// Returns Aerosol extinction optical depth for each layer and particle type: \f$ \tau_{aer\_ext,lp} \f$
    virtual ArrayAd<double, 2> aerosol_extinction_optical_depth_per_particle() const = 0;

    /// Returns Aerosol scattering optical depth for each layer and particle type: \f$ \tau_{aer\_sca,lp} \f$
    virtual ArrayAd<double, 2> aerosol_scattering_optical_depth_per_particle() const = 0;

    /// Returns aerosol phase function moments per particle
    /// Dimensions: num_moments x num_layers x num_scattering
    virtual const std::vector<ArrayAd<double, 3> > aerosol_phase_function_moments_per_particle(int num_moments = -1, int num_scattering = -1) const = 0;

    /// Gas Absorber optical depth summed over all particles:
    /// \f$ \tau_{gas,l} = \sum_{p=0}^{P} \tau_{gas,lp} \f$
    virtual ArrayAd<double, 1> gas_optical_depth_per_layer() const = 0;

    /// Aerosol extinction optical depth summed over all particles
    /// \f$ \tau_{aer\_ext,l} = \sum_{p=0}^{P} \tau_{aer\_ext,lp} \f$
    virtual ArrayAd<double, 1> aerosol_extinction_optical_depth_per_layer() const = 0;

    /// Aerosol scattering optical depth summed over all particles
    /// \f$ \tau_{aer\_sca,l} = \sum_{p=0}^{P} \tau_{aer\_sca,lp} \f$
    virtual ArrayAd<double, 1> aerosol_scattering_optical_depth_per_layer() const = 0;

    /// Total optical depth consisting of rayleigh + gas + aerosol
    /// Gas and aerosol contributions are allowed to possibly be zero.
    /// Computed as:
    /// \f$ \tau_{tot,l} = \tau_{ray,l} + \tau_{gas,l}  + \tau_{aer\_ext,l} \f$
    virtual ArrayAd<double, 1> total_optical_depth() const = 0;

    /// Total single scattering albedo
    /// Aerosol contribution is allowed to possibly be zero.
    /// Computed as:
    /// \f$ \omega_{tot,l} = \frac{ \tau_{ray,l} + \tau_{aer\_sca,l} }{ \tau_{tot,l} } \f$
    virtual ArrayAd<double, 1> total_single_scattering_albedo() const = 0;

    /// Fraction of the scattering attributable to rayleigh at each layer
    virtual ArrayAd<double, 1> rayleigh_fraction() const = 0;

    /// Fraction of the scattering attributable to each aerosol particle at each layer
    virtual ArrayAd<double, 2> aerosol_fraction() const = 0;

    /// Returns the portion of the total phase function attributable to
    /// aerosols. This is the summation over particles multiplied by the
    /// fraction of aerosol.
    ///
    /// Output dimensions: num_moments x num_layers x num_scattering
    virtual ArrayAd<double, 3> rayleigh_phase_function_moments_portion(int num_scattering) const = 0;

    /// Compute the portion of the total phase function attributable to
    /// rayleigh scattering. This returns the Rayleigh Greek Moments matrix
    /// multiplied by the fraction of rayleigh.
    ///
    /// Returned dimensions: num_moments x num_layers x num_scattering
    virtual ArrayAd<double, 3> aerosol_phase_function_moments_portion(int num_moments = -1, int num_scattering = -1) const = 0;

    /// Compute the total phase function moments which is the summation
    /// of the rayleigh and aerosol portions
    virtual ArrayAd<double, 3> total_phase_function_moments(int num_moments = -1, int num_scattering = -1) const = 0;

    /// Matrix that can be multiplied by the jacobians output by this class to provide jacobians with respect to the original input variables
    /// Size is num_layers x num_intermediate_jacobians x num_input_jacobians
    virtual blitz::Array<double, 3> intermediate_jacobian() const = 0;
    
};

}

#endif

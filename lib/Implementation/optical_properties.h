#ifndef PCA_OPTICAL_PROP_H
#define PCA_OPTICAL_PROP_H

#include "generic_object.h"

#include "array_ad.h"

#include "absorber.h"
#include "rayleigh.h"
#include "aerosol.h"

namespace FullPhysics {

/****************************************************************//**
  Represents the optical properties for a single spectral point.
  The intention is to contain all the information calculated from
  the atmosphere needed for a radiative transfer calculation along
  with intermediate computations useful by approximation methods.
 *******************************************************************/

class OpticalProperties : public virtual GenericObject {
public:

    OpticalProperties() : initialized(false) {}

    virtual void initialize(const ArrayAd<double, 1>& rayleigh_od, 
                            const ArrayAd<double, 2>& gas_od,
                            const ArrayAd<double, 2>& aerosol_ext_od,
                            const ArrayAd<double, 2>& aerosol_sca_od,
                            const std::vector<ArrayAd<double, 3> >& aerosol_pf_moments);

    virtual void initialize(const DoubleWithUnit spectral_point,
                            const int channel_index,
                            const boost::shared_ptr<Absorber>& absorber,
                            const boost::shared_ptr<Rayleigh>& rayleigh,
                            const boost::shared_ptr<Aerosol>& aerosol,
                            int num_pf_mom = -1, int num_scattering = -1);

    /// Deconstructor
    virtual ~OpticalProperties() = default;

    // Sizes of types of information stored

    virtual int number_layers() const { assert_sizes(); return rayleigh_optical_depth_.rows(); }
    virtual int number_gas_particles() const { assert_sizes(); return gas_optical_depth_per_particle_.cols(); }
    virtual int number_aerosol_particles() const { assert_sizes(); return aerosol_extinction_optical_depth_per_particle_.cols(); }
 
    // 
    // These accessors simply return what was passed in
    // 
    
    /// Returns Rayleigh optical depth for each layer: \f$ \tau_{ray,l} \f$
    virtual ArrayAd<double, 1> rayleigh_optical_depth() const { assert_init(); return rayleigh_optical_depth_; }

    /// Returns Gas Absorber optical depth for each layer and particle type: \f$ \tau_{gas,lp} \f$
    virtual ArrayAd<double, 2> gas_optical_depth_per_particle() const { assert_init(); return gas_optical_depth_per_particle_; }

    /// Returns Aerosol extinction optical depth for each layer and particle type: \f$ \tau_{aer\_ext,lp} \f$
    virtual ArrayAd<double, 2> aerosol_extinction_optical_depth_per_particle() const { assert_init(); return aerosol_extinction_optical_depth_per_particle_; }

    /// Returns Aerosol scattering optical depth for each layer and particle type: \f$ \tau_{aer\_sca,lp} \f$
    virtual ArrayAd<double, 2> aerosol_scattering_optical_depth_per_particle() const { assert_init(); return aerosol_scattering_optical_depth_per_particle_; }

    /// Returns aerosol phase function moments per particle
    /// Dimensions: num_moments x num_layers x num_scattering
    virtual const std::vector<ArrayAd<double, 3> >& aerosol_phase_function_moments_per_particle() const { assert_init(); return aerosol_phase_function_moments_per_particle_; }

    // 
    // These accessors only calculate their value if their stored value is empty
    // 

    virtual ArrayAd<double, 1> gas_optical_depth_per_layer() const;
    virtual ArrayAd<double, 1> aerosol_extinction_optical_depth_per_layer() const;
    virtual ArrayAd<double, 1> aerosol_scattering_optical_depth_per_layer() const;

    virtual ArrayAd<double, 1> total_optical_depth() const;
    virtual ArrayAd<double, 1> total_single_scattering_albedo() const;

    virtual ArrayAd<double, 1> rayleigh_fraction() const;
    virtual ArrayAd<double, 2> aerosol_fraction() const;

    virtual ArrayAd<double, 3> rayleigh_phase_function_moments_portion() const;
    virtual ArrayAd<double, 3> aerosol_phase_function_moments_portion() const;
    virtual ArrayAd<double, 3> total_phase_function_moments() const;

    /// Matrix that can be multiplied by the jacobians output by this class to provide jacobians with respect to the original input variables
    virtual blitz::Array<double, 3> intermediate_jacobian() const { assert_init(); return intermediate_jacobian_; }
    
protected:

    // Initialization protection
    bool initialized;
    void assert_init() const;
    void assert_sizes() const;

    virtual void initialize_with_jacobians(const ArrayAd<double, 1>& rayleigh_od, 
                                           const ArrayAd<double, 2>& gas_od,
                                           const ArrayAd<double, 2>& aerosol_ext_od,
                                           const ArrayAd<double, 2>& aerosol_sca_od,
                                           const std::vector<ArrayAd<double, 3> >& aerosol_pf_moments);

    // Dim: num_layers x num_gases
    ArrayAd<double, 2> gas_optical_depth_per_particle_;

    // Dim: num_layers
    ArrayAd<double, 1> rayleigh_optical_depth_;

    // Dim: num_layers x num_particles
    ArrayAd<double, 2> aerosol_extinction_optical_depth_per_particle_;
    ArrayAd<double, 2> aerosol_scattering_optical_depth_per_particle_;

    // Dim: num_moments x num_layers x num_scattering
    std::vector<ArrayAd<double, 3> > aerosol_phase_function_moments_per_particle_;

    // Dim: num_layers
    // Summation over all particles
    mutable ArrayAd<double, 1> gas_optical_depth_per_layer_;
    mutable ArrayAd<double, 1> aerosol_extinction_optical_depth_per_layer_;
    mutable ArrayAd<double, 1> aerosol_scattering_optical_depth_per_layer_;

    // Portion of total phase function attributable to aerosol
    mutable ArrayAd<double, 3> rayleigh_phase_function_moments_portion_;
    mutable ArrayAd<double, 3> aerosol_phase_function_moments_portion_;
    mutable ArrayAd<double, 3> total_phase_function_moments_;

    // Primary optical properties intended for radiative transfer
    // mutable since thise are computed on demand
    mutable ArrayAd<double, 1> total_optical_depth_;
    mutable ArrayAd<double, 1> total_single_scattering_albedo_;

    // Helper function and variable for a value used in total_single_scattering_albedo, rayleigh_fraction and aerosol_fraction
    ArrayAd<double, 1> scattering_sum() const;
    mutable ArrayAd<double, 1> scattering_sum_;

    // Convenient values needed in scattering moment calculation
    mutable ArrayAd<double, 1> rayleigh_fraction_;
    mutable ArrayAd<double, 2> aerosol_fraction_;
    
    // For conversion from jacobians wrt internal value to wrt input values
    blitz::Array<double, 3> intermediate_jacobian_;

    // Reflective surface
    ArrayAd<double, 1> surface_reflective_parameters_;

    // Thermal surface
    ArrayAd<double, 1> atmosphere_blackbody_;
    AutoDerivative<double> surface_blackbody_;

};

/****************************************************************//**
 *******************************************************************/

class OpticalPropertiesWrtRt : public virtual OpticalProperties {
public:

    OpticalPropertiesWrtRt() : OpticalProperties() {};

protected:

    virtual void initialize_with_jacobians(const ArrayAd<double, 1>& rayleigh_od, 
                                           const ArrayAd<double, 2>& gas_od,
                                           const ArrayAd<double, 2>& aerosol_ext_od,
                                           const ArrayAd<double, 2>& aerosol_sca_od,
                                           const std::vector<ArrayAd<double, 3> >& aerosol_pf_moments);

};

}

#endif

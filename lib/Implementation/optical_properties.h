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
  with intermediate computations useful in speed up methods.
 *******************************************************************/
class OpticalProperties : public virtual GenericObject {
public:

    OpticalProperties(const ArrayAd<double, 1>& rayleigh_od, 
                      const ArrayAd<double, 2>& gas_od,
                      const ArrayAd<double, 2>& aerosol_ext_od,
                      const ArrayAd<double, 2>& aerosol_sca_od);

    OpticalProperties(const DoubleWithUnit spectral_point,
                      const int channel_index,
                      const boost::shared_ptr<Absorber>& absorber,
                      const boost::shared_ptr<Rayleigh>& rayleigh,
                      const boost::shared_ptr<Aerosol>& aerosol);

    virtual ~OpticalProperties() = default;

    // These accessors simply return what was passed in
    virtual ArrayAd<double, 1> rayleigh_optical_depth() const { return rayleigh_optical_depth_; }
    virtual ArrayAd<double, 2> gas_optical_depth_per_particle() const { return gas_optical_depth_per_particle_; }
    virtual ArrayAd<double, 2> aerosol_extinction_optical_depth_per_particle() const { return aerosol_extinction_optical_depth_per_particle_; }
    virtual ArrayAd<double, 2> aerosol_scattering_optical_depth_per_particle() const { return aerosol_scattering_optical_depth_per_particle_; }

    // These accessors only calculate their value if their stored value is empty
    virtual ArrayAd<double, 1> gas_optical_depth_per_layer() const;
    virtual ArrayAd<double, 1> aerosol_extinction_optical_depth_per_layer() const;
    virtual ArrayAd<double, 1> aerosol_scattering_optical_depth_per_layer() const;

    virtual ArrayAd<double, 1> total_optical_depth() const;
    virtual ArrayAd<double, 1> total_single_scattering_albedo() const;

private:

    // Intermediate optical properties used in the computation of primary properties

    // Dim: num_layers x num_gases
    ArrayAd<double, 2> gas_optical_depth_per_particle_;

    // Dim: num_layers
    ArrayAd<double, 1> rayleigh_optical_depth_;

    // Dim: num_layers x num_particles
    ArrayAd<double, 2> aerosol_extinction_optical_depth_per_particle_;
    ArrayAd<double, 2> aerosol_scattering_optical_depth_per_particle_;

    // Dim: num_layers
    // Summation over all particles
    mutable ArrayAd<double, 1> gas_optical_depth_per_layer_;
    mutable ArrayAd<double, 1> aerosol_extinction_optical_depth_per_layer_;
    mutable ArrayAd<double, 1> aerosol_scattering_optical_depth_per_layer_;

    // Primary optical properties intended for radiative transfer
    // mutable since thise are computed on demand
    mutable ArrayAd<double, 1> total_optical_depth_;
    mutable ArrayAd<double, 1> total_single_scattering_albedo_;

    // Reflective surface
    ArrayAd<double, 1> surface_reflective_parameters_;

    // Thermal surface
    ArrayAd<double, 1> atmosphere_blackbody_;
    AutoDerivative<double> surface_blackbody_;

    // Convenient values needed in scattering moment calculation
    ArrayAd<double, 1> rayleigh_fraction_;
    ArrayAd<double, 2> aerosol_fraction_;


};

}

#endif

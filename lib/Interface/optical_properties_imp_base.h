#ifndef OPTICAL_PROP_IMP_BASE_H
#define OPTICAL_PROP_IMP_BASE_H

#include "optical_properties.h"

#include "double_with_unit.h"
#include "array_ad.h"

#include "aerosol.h"

namespace FullPhysics {

/****************************************************************//**
  Helper class to return aerosol phase functions. Due to the 
  computational expense of computing these values it is a good idea
  to delay computation until the number of moments and scattering
  matrix elements is known. These classes provide wrappers for 
  encapsulating those calls.
 *******************************************************************/

class AerosolPhaseFunctionHelper :
    public Printable<AerosolPhaseFunctionHelper> {
public:
  AerosolPhaseFunctionHelper() {}
  virtual ~AerosolPhaseFunctionHelper() {}
  virtual const std::vector<ArrayAd<double, 3> >
  phase_function_moments_per_particle(int num_moments = -1,
				      int num_scattering = -1) const = 0;
//-----------------------------------------------------------------------
/// Print description of object.
//-----------------------------------------------------------------------

  virtual void print(std::ostream& Os) const
  {Os << "AerosolPhaseFunctionHelper";}
private:
  friend class boost::serialization::access;
  template<class Archive>
  void serialize(Archive & ar, const unsigned int version);
};

/****************************************************************//**
  Pass through passed aerosol moments, in this case the number of
  moments and scattering are trimmed from the inputs
*******************************************************************/

class AerosolPhaseFunctionPassThruHelper : public AerosolPhaseFunctionHelper {
public:
  AerosolPhaseFunctionPassThruHelper
  (const std::vector<ArrayAd<double, 3> >& aerosol_pf_moments);
  virtual const std::vector<ArrayAd<double, 3> >
  phase_function_moments_per_particle(int num_moments = -1,
				      int num_scattering = -1) const;
  virtual void print(std::ostream& Os) const
  {Os << "AerosolPhaseFunctionPassThruHelper";}
private:
  std::vector<ArrayAd<double, 3> > aerosol_pf_moments_in;
  AerosolPhaseFunctionPassThruHelper() {}
  friend class boost::serialization::access;
  template<class Archive>
  void serialize(Archive & ar, const unsigned int version);
};

/****************************************************************//**
  Compute aerosol phase functions using the Aerosol class hiearchy
*******************************************************************/

class AerosolPhaseFunctionComputeHelper : public AerosolPhaseFunctionHelper {
public:
  AerosolPhaseFunctionComputeHelper(const DoubleWithUnit spectral_point,
				    const boost::shared_ptr<Aerosol>& aerosol);
  virtual const std::vector<ArrayAd<double, 3> >
  phase_function_moments_per_particle(int num_moments = -1,
				      int num_scattering = -1) const;
  virtual void print(std::ostream& Os) const
  {Os << "AerosolPhaseFunctionComputeHelper";}
private:
  double wn;
  boost::shared_ptr<Aerosol> aerosol_;
  AerosolPhaseFunctionComputeHelper() {}
  friend class boost::serialization::access;
  template<class Archive>
  void serialize(Archive & ar, const unsigned int version);
};

/****************************************************************//**
  Represents the optical properties for a single spectral point.
  The intention is to contain all the information calculated from
  the atmosphere needed for a radiative transfer calculation along
  with intermediate computations useful by approximation methods.
 *******************************************************************/

class OpticalPropertiesImpBase : public virtual OpticalProperties {
public:
  OpticalPropertiesImpBase() : initialized(false), cached_num_moments(-1),
			       cached_num_scattering(-1) {}

//-----------------------------------------------------------------------
/// Deconstructor
//-----------------------------------------------------------------------

  virtual ~OpticalPropertiesImpBase() = default;

//-----------------------------------------------------------------------
/// Sizes of types of information stored
//-----------------------------------------------------------------------

  virtual int number_layers() const
  { assert_sizes(); return rayleigh_optical_depth_.rows(); }
  virtual int number_gas_particles() const
  { assert_sizes(); return gas_optical_depth_per_particle_.cols(); }
  virtual int number_aerosol_particles() const
  {
    assert_sizes();
    return aerosol_extinction_optical_depth_per_particle_.cols();
  }
 
//-----------------------------------------------------------------------
// These accessors simply return what was passed in
//-----------------------------------------------------------------------
    
//-----------------------------------------------------------------------
/// Returns Rayleigh optical depth for each layer: \f$ \tau_{ray,l}
/// \f$
//-----------------------------------------------------------------------
  
  virtual ArrayAd<double, 1> rayleigh_optical_depth() const
  { assert_init(); return rayleigh_optical_depth_; }

//-----------------------------------------------------------------------
/// Returns Gas Absorber optical depth for each layer and particle
/// type: \f$ \tau_{gas,lp} \f$
//-----------------------------------------------------------------------
  
  virtual ArrayAd<double, 2> gas_optical_depth_per_particle() const
  { assert_init(); return gas_optical_depth_per_particle_; }

//-----------------------------------------------------------------------
/// Returns Aerosol extinction optical depth for each layer and particle type: \f$ \tau_{aer\_ext,lp} \f$
//-----------------------------------------------------------------------

  virtual ArrayAd<double, 2> aerosol_extinction_optical_depth_per_particle()
    const
  { assert_init(); return aerosol_extinction_optical_depth_per_particle_; }

//-----------------------------------------------------------------------
/// Returns Aerosol scattering optical depth for each layer and particle type: \f$ \tau_{aer\_sca,lp} \f$
//-----------------------------------------------------------------------

  virtual ArrayAd<double, 2> aerosol_scattering_optical_depth_per_particle()
    const
  { assert_init(); return aerosol_scattering_optical_depth_per_particle_; }

//-----------------------------------------------------------------------
/// Returns aerosol phase function moments per particle
/// Dimensions: num_moments x num_layers x num_scattering
//-----------------------------------------------------------------------
  
  virtual const std::vector<ArrayAd<double, 3> >
  aerosol_phase_function_moments_per_particle(int num_moments = -1,
					      int num_scattering = -1) const;

//-----------------------------------------------------------------------
// These accessors only calculate their value if their stored value is empty
//-----------------------------------------------------------------------

  virtual ArrayAd<double, 1> gas_optical_depth_per_layer() const;
  virtual ArrayAd<double, 1> aerosol_extinction_optical_depth_per_layer() const;
  virtual ArrayAd<double, 1> aerosol_scattering_optical_depth_per_layer() const;

  virtual ArrayAd<double, 1> total_optical_depth() const;
  virtual ArrayAd<double, 1> total_single_scattering_albedo() const;

  virtual ArrayAd<double, 1> rayleigh_fraction() const;
  virtual ArrayAd<double, 2> aerosol_fraction() const;

  virtual boost::shared_ptr<AerosolPhaseFunctionHelper>
  aerosol_phase_function_helper() const
  { return aerosol_phase_function_helper_; }

  virtual ArrayAd<double, 3>
  rayleigh_phase_function_moments_portion(int num_scattering = -1) const;
  virtual ArrayAd<double, 3>
  aerosol_phase_function_moments_portion(int num_moments = -1,
					 int num_scattering = -1) const;
  virtual ArrayAd<double, 3>
  total_phase_function_moments(int num_moments = -1,
			       int num_scattering = -1) const;

//-----------------------------------------------------------------------
/// Matrix that can be multiplied by the jacobians output by this
/// class to provide jacobians with respect to the original input
/// variables Size is num_layers x num_intermediate_jacobians x
/// num_input_jacobians
//-----------------------------------------------------------------------

  virtual blitz::Array<double, 3> intermediate_jacobian() const
  { assert_init(); return intermediate_jacobian_; }
    
protected:
  // Initialization protection
  bool initialized;
  void assert_init() const;
  void assert_sizes() const;

  // Dim: num_layers x num_gases
  ArrayAd<double, 2> gas_optical_depth_per_particle_;

  // Dim: num_layers
  ArrayAd<double, 1> rayleigh_optical_depth_;

  // Dim: num_layers x num_particles
  ArrayAd<double, 2> aerosol_extinction_optical_depth_per_particle_;
  ArrayAd<double, 2> aerosol_scattering_optical_depth_per_particle_;

  // Dim: num_moments x num_layers x num_scattering
  mutable int cached_num_moments;
  mutable int cached_num_scattering;
  boost::shared_ptr<AerosolPhaseFunctionHelper> aerosol_phase_function_helper_;
  mutable std::vector<ArrayAd<double, 3> >
  aerosol_phase_function_moments_per_particle_;

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

  // Helper function and variable for a value used in
  // total_single_scattering_albedo, rayleigh_fraction and
  // aerosol_fraction
  ArrayAd<double, 1> scattering_sum() const;
  mutable ArrayAd<double, 1> scattering_sum_;

  // Convenient values needed in scattering moment calculation
  mutable ArrayAd<double, 1> rayleigh_fraction_;
  mutable ArrayAd<double, 2> aerosol_fraction_;
    
  // For conversion from jacobians wrt internal value to wrt input values
  // Size is num_layers x num_intermediate_jacobians x num_input_jacobians
  blitz::Array<double, 3> intermediate_jacobian_;

  // Reflective surface
  ArrayAd<double, 1> surface_reflective_parameters_;

  // Thermal surface
  ArrayAd<double, 1> atmosphere_blackbody_;
  AutoDerivative<double> surface_blackbody_;
private:  
  friend class boost::serialization::access;
  template<class Archive>
  void serialize(Archive & ar, const unsigned int version);
};

}

FP_EXPORT_KEY(AerosolPhaseFunctionHelper);
FP_EXPORT_KEY(AerosolPhaseFunctionPassThruHelper);
FP_EXPORT_KEY(AerosolPhaseFunctionComputeHelper);
FP_EXPORT_KEY(OpticalPropertiesImpBase);
#endif

// -*- mode: c++; -*-
// (Not really c++, but closest emacs mode)

%include "fp_common.i"

%{
#include "optical_properties_imp_base.h"
#include "sub_state_vector_array.h"
#include "rayleigh.h"
#include "altitude.h"
%}

%base_import(optical_properties)

%import "array_ad.i"
%import "double_with_unit.i"

%import "aerosol.i"

%fp_shared_ptr(FullPhysics::AerosolPhaseFunctionHelper)
%fp_shared_ptr(FullPhysics::AerosolPhaseFunctionPassThruHelper)
%fp_shared_ptr(FullPhysics::AerosolPhaseFunctionComputeHelper)
%fp_shared_ptr(FullPhysics::OpticalPropertiesImpBase)

namespace FullPhysics {
class AerosolPhaseFunctionHelper : public GenericObject {
public:
  virtual const std::vector<ArrayAd<double, 3> >
  phase_function_moments_per_particle(int num_moments = -1,
				      int num_scattering = -1) const = 0;
  std::string print_to_string() const;
  std::string print_parent() const;
  %pickle_serialization();
};

class AerosolPhaseFunctionPassThruHelper : public AerosolPhaseFunctionHelper {
public:
  AerosolPhaseFunctionPassThruHelper
  (const std::vector<ArrayAd<double, 3> >& aerosol_pf_moments);
  virtual const std::vector<ArrayAd<double, 3> >
  phase_function_moments_per_particle(int num_moments = -1,
				      int num_scattering = -1) const;
  %pickle_serialization();
};

class AerosolPhaseFunctionComputeHelper : public AerosolPhaseFunctionHelper {
public:
  AerosolPhaseFunctionComputeHelper(const DoubleWithUnit spectral_point,
				    const boost::shared_ptr<Aerosol>& aerosol);
  virtual const std::vector<ArrayAd<double, 3> >
  phase_function_moments_per_particle(int num_moments = -1,
				      int num_scattering = -1) const;
  %pickle_serialization();
};

class OpticalPropertiesImpBase : public OpticalProperties {
public:
  OpticalPropertiesImpBase();
  virtual int number_layers() const;
  virtual int number_gas_particles() const;
  virtual int number_aerosol_particles() const;
  virtual ArrayAd<double, 1> rayleigh_optical_depth() const;
  virtual ArrayAd<double, 2> gas_optical_depth_per_particle() const;
  virtual ArrayAd<double, 2> aerosol_extinction_optical_depth_per_particle()
    const;
  virtual ArrayAd<double, 2> aerosol_scattering_optical_depth_per_particle()
    const;
  virtual const std::vector<ArrayAd<double, 3> >
  aerosol_phase_function_moments_per_particle(int num_moments = -1,
					      int num_scattering = -1) const;
  virtual ArrayAd<double, 1> gas_optical_depth_per_layer() const;
  virtual ArrayAd<double, 1> aerosol_extinction_optical_depth_per_layer() const;
  virtual ArrayAd<double, 1> aerosol_scattering_optical_depth_per_layer() const;
  virtual ArrayAd<double, 1> total_optical_depth() const;
  virtual ArrayAd<double, 1> total_single_scattering_albedo() const;
  virtual ArrayAd<double, 1> rayleigh_fraction() const;
  virtual ArrayAd<double, 2> aerosol_fraction() const;
  virtual boost::shared_ptr<AerosolPhaseFunctionHelper>
  aerosol_phase_function_helper() const;
  virtual ArrayAd<double, 3>
  rayleigh_phase_function_moments_portion(int num_scattering = -1) const;
  virtual ArrayAd<double, 3>
  aerosol_phase_function_moments_portion(int num_moments = -1,
					 int num_scattering = -1) const;
  virtual ArrayAd<double, 3>
  total_phase_function_moments(int num_moments = -1,
			       int num_scattering = -1) const;
  virtual blitz::Array<double, 3> intermediate_jacobian() const;
  %pickle_serialization();
};
  
}  

